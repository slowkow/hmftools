//! The files Virus Interpreter reads and writes: the taxonomy, reporting and blacklisting databases
//! (TaxonomyDbFile, VirusReportingDbFile, VirusBlacklistingDbFile), VIRUSBreakend's summary
//! (hmf-common VirusBreakendFile) and the annotated viruses (hmf-common AnnotatedVirusFile).

use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, bail, Context, Result};
use hmf_common::fmt::java_double_to_string;
use hmf_common::java_text::{parse_double, parse_int, split};

use hmf_common::warn;

/// String.trim(): characters up to ' ' removed from both ends
pub fn java_trim(s: &str) -> &str {
    s.trim_matches(|c: char| c <= ' ')
}

fn int(v: &str) -> Result<i32> {
    parse_int(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
}

fn dbl(v: &str) -> Result<f64> {
    parse_double(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
}

/// Files.readAllLines
fn read_all_lines(path: &str) -> Result<Vec<String>> {
    hmf_common::java_text::read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))
}

// ---- VirusType, VirusLikelihoodType, VirusBreakendQCStatus ----

/// VirusType: its name as written ("HHV-8" for HHV8)
pub const VIRUS_TYPES: [&str; 5] = ["MCV", "EBV", "HPV", "HBV", "HHV-8"];

/// VirusType.fromVirusName
pub fn virus_type_from_name(name: &str) -> Result<&'static str> {
    VIRUS_TYPES
        .iter()
        .find(|t| **t == name)
        .copied()
        .ok_or_else(|| anyhow!("java.lang.IllegalStateException: Cannot resolve virus name: {name}"))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VirusLikelihoodType {
    High,
    Low,
    Unknown,
}

impl VirusLikelihoodType {
    /// valueOf(name)
    pub fn value_of(s: &str) -> Result<VirusLikelihoodType> {
        Ok(match s {
            "HIGH" => VirusLikelihoodType::High,
            "LOW" => VirusLikelihoodType::Low,
            "UNKNOWN" => VirusLikelihoodType::Unknown,
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.virusinterpreter.VirusLikelihoodType.{s}"),
        })
    }
    /// the enum name, as String.valueOf writes it
    pub fn name(self) -> &'static str {
        match self {
            VirusLikelihoodType::High => "HIGH",
            VirusLikelihoodType::Low => "LOW",
            VirusLikelihoodType::Unknown => "UNKNOWN",
        }
    }
}

/// VirusBreakendQCStatus, in enum order
pub const QC_STATUSES: [&str; 6] =
    ["NO_ABNORMALITIES", "LOW_VIRAL_COVERAGE", "EXCESSIVE_VIRAL_COVERAGE", "ASSEMBLY_DOWNSAMPLED", "CHILD_TAXID_REFERENCE", "UNCLEAR_TAXID_ASSIGNMENT"];

/// VirusBreakendQCStatus.convert: NO_ABNORMALITIES for an empty value
pub fn qc_status_convert(s: &str) -> Result<&'static str> {
    if s.is_empty() {
        return Ok(QC_STATUSES[0]);
    }
    QC_STATUSES
        .iter()
        .find(|q| **q == s)
        .copied()
        .ok_or_else(|| anyhow!("java.lang.IllegalStateException: Cannot resolve QC Status of VIRUSBreakend: {s}"))
}

// ---- the taxonomy database ----

/// TaxonomyDb, holding only the taxids the run looks up
#[derive(Clone, Debug, Default)]
pub struct TaxonomyDb {
    pub names: HashMap<i32, String>,
    /// the number of taxids in the file (TaxonomyDb.count)
    pub count: usize,
}

impl TaxonomyDb {
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn taxid_exists(&self, taxid: i32) -> bool {
        self.names.contains_key(&taxid)
    }

    /// lookupName: "" (and a warning) for a taxid not in the database
    pub fn lookup_name(&self, taxid: i32) -> String {
        match self.names.get(&taxid) {
            Some(n) => n.clone(),
            None => {
                warn!("Could not match taxid '{}' to a name in taxonomy DB", taxid);
                String::new()
            }
        }
    }
}

/// TaxonomyDbFile.loadFromTsv: every line (the file has no header) of two tab-separated fields (trailing
/// empty fields dropped, as String.split does) is a taxid and its name, both trimmed, a later line replacing
/// an earlier one; other lines are warned about. A taxid that is not a number ends the run, and so does text
/// that is not UTF-8 (Files.readAllLines), as in Java.
///
/// Java keeps the whole file (2.6 million names) in a HashMap; the port reads it in blocks and keeps the
/// names of `wanted` only (all of them if None, when `count` is the number of taxids as TaxonomyDb.count).
pub fn load_taxonomy_db(path: &str, wanted: Option<&HashSet<i32>>) -> Result<TaxonomyDb> {
    use std::io::BufRead as _;
    let file = std::fs::File::open(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let mut reader = std::io::BufReader::with_capacity(4 << 20, file);
    let mut db = TaxonomyDb::default();
    let mut all: HashSet<i32> = HashSet::new();
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    loop {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf).with_context(|| format!("reading {path}"))?;
        if n == 0 {
            break;
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
        let text = std::str::from_utf8(&buf).map_err(|_| anyhow!("java.nio.charset.MalformedInputException: Input length = 1 ({path})"))?;
        // BufferedReader also ends a line at a lone '\r'
        for line in text.split('\r') {
            let trimmed = line.trim_end_matches('\t');
            let mut tabs = memchr::memchr_iter(b'\t', trimmed.as_bytes());
            let first = tabs.next();
            if let (Some(t), None) = (first, tabs.next()) {
                let taxid = int(java_trim(&trimmed[..t]))?;
                if wanted.is_none() {
                    all.insert(taxid);
                }
                if wanted.is_none_or(|w| w.contains(&taxid)) {
                    db.names.insert(taxid, java_trim(&trimmed[t + 1..]).to_string());
                }
            } else {
                warn!("Suspicious line detected in taxonomy db tsv: {}", line);
            }
        }
    }
    db.count = all.len();
    Ok(db)
}

// ---- the blacklisting and reporting databases ----

/// VirusBlacklistingDbFile.loadFromTsv: the first field of each line of four fields, after the header
pub fn load_blacklisting_db(path: &str) -> Result<Vec<i32>> {
    let lines = read_all_lines(path)?;
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let parts = split(line, '\t');
        if parts.len() == 4 {
            out.push(int(java_trim(parts[0]))?);
        } else {
            warn!("Suspicious line detected in virus blacklisting db tsv: {}", line);
        }
    }
    Ok(out)
}

/// VirusReportingDb
#[derive(Clone, Debug, PartialEq)]
pub struct VirusReportingDb {
    pub virus_interpretation: &'static str,
    pub integrated_minimal_coverage: Option<i32>,
    pub non_integrated_minimal_coverage: Option<i32>,
    pub virus_driver_likelihood_type: VirusLikelihoodType,
}

/// VirusReportingDbModel: species taxid to its reporting
#[derive(Clone, Debug, Default)]
pub struct VirusReportingDbModel {
    pub species: HashMap<i32, VirusReportingDb>,
}

impl VirusReportingDbModel {
    pub fn has_interpretation(&self, species_taxid: i32) -> bool {
        self.species.contains_key(&species_taxid)
    }
    fn note_missing(&self, species_taxid: i32) {
        if !self.has_interpretation(species_taxid) {
            hmf_common::debug!("No interpretation found for virus with species taxid {}", species_taxid);
        }
    }
    pub fn interpret_virus_species(&self, species_taxid: i32) -> Option<&'static str> {
        self.note_missing(species_taxid);
        self.species.get(&species_taxid).map(|r| r.virus_interpretation)
    }
    pub fn virus_likelihood_type(&self, species_taxid: i32) -> VirusLikelihoodType {
        self.note_missing(species_taxid);
        self.species.get(&species_taxid).map_or(VirusLikelihoodType::Unknown, |r| r.virus_driver_likelihood_type)
    }
    pub fn integrated_minimal_coverage(&self, species_taxid: i32) -> Option<i32> {
        self.note_missing(species_taxid);
        self.species.get(&species_taxid).and_then(|r| r.integrated_minimal_coverage)
    }
    pub fn non_integrated_minimal_coverage(&self, species_taxid: i32) -> Option<i32> {
        self.note_missing(species_taxid);
        self.species.get(&species_taxid).and_then(|r| r.non_integrated_minimal_coverage)
    }
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn count(&self) -> usize {
        self.species.len()
    }
}

/// VirusReportingDbFile.buildFromTsv: lines of six fields after the header (species taxid, virus name, the
/// integrated and non-integrated minimal coverages, empty for none, the driver likelihood and a remark); a
/// later line of a species replaces an earlier one
pub fn load_reporting_db(path: &str) -> Result<VirusReportingDbModel> {
    let lines = read_all_lines(path)?;
    let mut model = VirusReportingDbModel::default();
    for line in lines.iter().skip(1) {
        let parts = split(line, '\t');
        if parts.len() == 6 {
            let species_taxid = int(java_trim(parts[0]))?;
            let coverage = |s: &str| -> Result<Option<i32>> {
                let t = java_trim(s);
                if t.is_empty() { Ok(None) } else { int(t).map(Some) }
            };
            let virus_interpretation = virus_type_from_name(java_trim(parts[1]))?;
            let integrated_minimal_coverage = coverage(parts[2])?;
            let non_integrated_minimal_coverage = coverage(parts[3])?;
            let virus_driver_likelihood_type = VirusLikelihoodType::value_of(parts[4])?;
            model.species.insert(
                species_taxid,
                VirusReportingDb { virus_interpretation, integrated_minimal_coverage, non_integrated_minimal_coverage, virus_driver_likelihood_type },
            );
        } else {
            warn!("Suspicious line detected in virus reporting db tsv: {}", line);
        }
    }
    Ok(model)
}

// ---- VIRUSBreakend's summary ----

/// VirusBreakend: the fields Virus Interpreter reads (every field is parsed, as Java does)
#[derive(Clone, Debug, PartialEq)]
pub struct VirusBreakend {
    pub taxid_genus: i32,
    pub taxid_species: i32,
    pub reference_taxid: i32,
    pub coverage: f64,
    pub mean_depth: f64,
    pub integrations: i32,
    pub qc_status: &'static str,
}

impl Default for VirusBreakend {
    fn default() -> Self {
        VirusBreakend { taxid_genus: 0, taxid_species: 0, reference_taxid: 0, coverage: 0.0, mean_depth: 0.0, integrations: 0, qc_status: QC_STATUSES[0] }
    }
}

/// VirusBreakendFile.read: each line after the first; the fields of a line split at tabs with trailing empty
/// fields dropped (so an empty QC status, the last field, reads as NO_ABNORMALITIES)
pub fn read_virus_breakends(path: &str) -> Result<Vec<VirusBreakend>> {
    let lines = read_all_lines(path)?;
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let v = split(line, '\t');
        let get = |i: usize| -> Result<&str> {
            v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))
        };
        // the fields in Java's order of parsing, so that the first bad one is the one reported
        let taxid_genus = int(get(0)?)?;
        get(1)?;
        int(get(2)?)?;
        let taxid_species = int(get(3)?)?;
        get(4)?;
        int(get(5)?)?;
        int(get(6)?)?;
        get(7)?;
        int(get(8)?)?;
        int(get(9)?)?;
        get(10)?;
        let reference_taxid = int(get(11)?)?;
        for i in 12..14 {
            int(get(i)?)?;
        }
        get(14)?;
        for i in 15..19 {
            int(get(i)?)?;
        }
        let coverage = dbl(get(19)?)?;
        let mean_depth = dbl(get(20)?)?;
        dbl(get(21)?)?;
        dbl(get(22)?)?;
        let integrations = int(get(23)?)?;
        let qc_status = if v.len() == 25 { qc_status_convert(v[24])? } else { QC_STATUSES[0] };
        out.push(VirusBreakend { taxid_genus, taxid_species, reference_taxid, coverage, mean_depth, integrations, qc_status });
    }
    Ok(out)
}

// ---- the annotated viruses ----

/// AnnotatedVirus
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotatedVirus {
    pub taxid: i32,
    pub name: String,
    pub qc_status: &'static str,
    pub integrations: i32,
    pub interpretation: Option<&'static str>,
    pub percentage_covered: f64,
    pub mean_coverage: f64,
    pub expected_clonal_coverage: Option<f64>,
    pub reported: bool,
    pub blacklisted: Option<bool>,
    pub virus_driver_likelihood_type: VirusLikelihoodType,
}

/// AnnotatedVirusFile.generateFileName
pub fn annotated_virus_filename(output_prefix: &str, sample_id: &str) -> String {
    format!("{output_prefix}{sample_id}.virus.annotated.tsv")
}

pub const ANNOTATED_VIRUS_HEADER: &str =
    "taxid\tname\tqcStatus\tintegrations\tinterpretation\tpercentageCovered\tmeanCoverage\texpectedClonalCoverage\treported\tblacklisted\tdriverLikelihood";

/// AnnotatedVirusFile.toLines, each line ended by a newline (Files.write)
pub fn annotated_virus_text(viruses: &[AnnotatedVirus]) -> String {
    let mut out = String::from(ANNOTATED_VIRUS_HEADER);
    out.push('\n');
    for v in viruses {
        let fields = [
            v.taxid.to_string(),
            v.name.clone(),
            v.qc_status.to_string(),
            v.integrations.to_string(),
            v.interpretation.unwrap_or("null").to_string(),
            java_double_to_string(v.percentage_covered),
            java_double_to_string(v.mean_coverage),
            v.expected_clonal_coverage.map_or("null".to_string(), java_double_to_string),
            v.reported.to_string(),
            v.blacklisted.map_or("null".to_string(), |b| b.to_string()),
            v.virus_driver_likelihood_type.name().to_string(),
        ];
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}
