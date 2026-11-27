//! Readers of the result files of the smaller tools, as hmf-common reads them (ORANGE): CHORD's
//! prediction (ChordDataFile), PEACH's best haplotypes (PeachGenotypeFile, with the old Python format),
//! SIGS' allocations (SignatureAllocationFile) and the signature etiologies (SnvSigUtils, the resource
//! `sigs/signatures_etiology.tsv`), Virus Interpreter's annotated viruses (AnnotatedVirusFile), CUPPA's
//! visualisation data (CuppaPredictions.fromTsv), and the pipeline version file (PipelineVersionFile).
//! Errors where Java throws are `Err`.

use anyhow::{anyhow, bail, Context, Result};

use crate::java_hash::JavaHashMap;
use crate::java_text::{fields_index_map, read_lines, split, split_all};
use crate::linx_files::{enum_value, parse_double_value, parse_int_value, Row};

fn dir(base_path: &str) -> String {
    if base_path.ends_with('/') { base_path.to_string() } else { format!("{base_path}/") }
}

fn lines_of(path: &str) -> Result<Vec<String>> {
    read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))
}

fn header(lines: &[String], path: &str) -> Result<std::collections::HashMap<String, usize>> {
    let h = lines.first().ok_or_else(|| anyhow!("{path}: java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    Ok(fields_index_map(h, '\t'))
}

// ---------------------------------------------------------------------------------------------
// CHORD

#[derive(Clone, Debug, PartialEq)]
pub struct ChordData {
    pub brca1_value: f64,
    pub brca2_value: f64,
    pub hrd_value: f64,
    /// ChordStatus name: CANNOT_BE_DETERMINED, HR_PROFICIENT, HR_DEFICIENT or UNKNOWN
    pub hr_status: String,
    pub hrd_type: String,
    pub remarks_hr_status: String,
    pub remarks_hrd_type: String,
}

/// ChordDataFile.generateFilename
pub fn chord_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.chord.prediction.tsv", dir(base_path))
}

/// ChordDataFile.extractHrStatus (an unknown value is UNKNOWN, with a warning)
pub fn chord_hr_status(s: &str) -> String {
    match s {
        "cannot_be_determined" => "CANNOT_BE_DETERMINED",
        "HR_proficient" => "HR_PROFICIENT",
        "HR_deficient" => "HR_DEFICIENT",
        _ => {
            crate::warn!("Unknown CHORD HR status: '{s}'");
            "UNKNOWN"
        }
    }
    .to_string()
}

/// ChordDataFile.read: exactly a header and one data line
pub fn read_chord(path: &str) -> Result<ChordData> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    if lines.len() != 2 {
        bail!("invalid chord file({path}), expected 2 lines containing header and data");
    }
    let f = fields_index_map(&lines[0], '\t');
    let values = split_all(&lines[1], '\t');
    let n = values.len();
    let r = Row { file: path, fields: &f, values };
    let remarks = |name: &str| -> Result<String> {
        Ok(if n >= 8 && f.contains_key(name) { r.string(name)? } else { String::new() })
    };
    let remarks_hr_status = remarks("remarks_hr_status")?;
    let remarks_hrd_type = remarks("remarks_hrd_type")?;
    Ok(ChordData {
        brca1_value: r.double("p_BRCA1")?,
        brca2_value: r.double("p_BRCA2")?,
        hrd_value: r.double("p_hrd")?,
        hr_status: chord_hr_status(r.str("hr_status")?),
        hrd_type: r.string("hrd_type")?,
        remarks_hr_status,
        remarks_hrd_type,
    })
}

// ---------------------------------------------------------------------------------------------
// PEACH

#[derive(Clone, Debug, PartialEq)]
pub struct PeachGenotype {
    pub gene: String,
    pub allele: String,
    pub allele_count: i32,
    pub function: String,
    pub linked_drugs: String,
    pub url_prescription_info: String,
    pub panel_version: String,
    pub repo_version: String,
}

impl PeachGenotype {
    /// PeachGenotype.haplotype(): PeachUtil.convertToHaplotypeString (Err for a count other than 1 or 2,
    /// as Java throws when the value is derived)
    pub fn haplotype(&self) -> Result<String> {
        if self.allele == "Unresolved Haplotype" {
            return Ok(self.allele.clone());
        }
        let zygosity = match self.allele_count {
            1 => "HET",
            2 => "HOM",
            n => bail!("java.lang.IllegalArgumentException: Could not convert allele count {n} to a zygosity"),
        };
        Ok(format!("{}_{zygosity}", self.allele))
    }
}

/// PeachGenotypeFile.generateFileName
pub fn peach_genotype_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.peach.haplotypes.best.tsv", dir(base_path))
}

const PEACH_PYTHON_HEADER: &str =
    "gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\thaplotype_only\tzygosity_only";

/// PeachGenotypeFile.read: the Java PEACH file, or the old Python PEACH file (recognised by its header)
pub fn read_peach_genotypes(path: &str) -> Result<Vec<PeachGenotype>> {
    let lines = lines_of(path)?;
    let f = header(&lines, path)?;
    let python = lines[0] == PEACH_PYTHON_HEADER;
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        if python {
            let zygosity = r.str("zygosity_only")?;
            let allele_count = match zygosity {
                "HET" => 1,
                "HOM" | "N/A" => 2,
                z => bail!("java.lang.IllegalArgumentException: Could not convert zygosity {z} to an allele count"),
            };
            out.push(PeachGenotype {
                gene: r.string("gene")?,
                allele: r.string("haplotype_only")?,
                allele_count,
                function: r.string("function")?,
                linked_drugs: r.string("linked_drugs")?,
                url_prescription_info: r.string("url_prescription_info")?,
                panel_version: r.string("panel_version")?,
                repo_version: r.string("repo_version")?,
            });
        } else {
            out.push(PeachGenotype {
                gene: r.string("gene")?,
                allele: r.string("haplotype")?,
                allele_count: r.int("count")?,
                function: r.string("function")?,
                linked_drugs: r.string("linkedDrugs")?,
                url_prescription_info: r.string("prescriptionUrls")?,
                panel_version: String::new(),
                repo_version: String::new(),
            });
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// SIGS

#[derive(Clone, Debug, PartialEq)]
pub struct SignatureAllocation {
    pub signature: String,
    pub allocation: f64,
    pub percent: f64,
}

/// SignatureAllocationFile.generateFilename
pub fn signature_allocation_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.sig.allocation.tsv", dir(base_path))
}

/// SignatureAllocationFile.read
pub fn read_signature_allocations(path: &str) -> Result<Vec<SignatureAllocation>> {
    let lines = lines_of(path)?;
    let f = header(&lines, path)?;
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        out.push(SignatureAllocation { signature: r.string("signature")?, allocation: r.double("allocation")?, percent: r.double("percent")? });
    }
    Ok(out)
}

const SIGNATURE_ETIOLOGIES: &str = include_str!("../resources/sigs/signatures_etiology.tsv");

/// SnvSigUtils.loadSnvSignatureEtiologies: signature name to etiology (a HashMap in Java)
pub fn snv_signature_etiologies() -> JavaHashMap<String, String> {
    let lines: Vec<&str> = SIGNATURE_ETIOLOGIES.lines().collect();
    let f = fields_index_map(lines[0], '\t');
    let (s, e) = (f["signature"], f["etiology"]);
    let mut m = JavaHashMap::new();
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        m.insert(v[s].to_string(), v[e].to_string());
    }
    m
}

// ---------------------------------------------------------------------------------------------
// Virus Interpreter

/// VirusBreakendQCStatus names
pub const VIRUS_BREAKEND_QC_STATUSES: [&str; 6] = [
    "NO_ABNORMALITIES",
    "LOW_VIRAL_COVERAGE",
    "EXCESSIVE_VIRAL_COVERAGE",
    "ASSEMBLY_DOWNSAMPLED",
    "CHILD_TAXID_REFERENCE",
    "UNCLEAR_TAXID_ASSIGNMENT",
];

#[derive(Clone, Debug, PartialEq)]
pub struct AnnotatedVirus {
    pub taxid: i32,
    pub name: String,
    pub qc_status: String,
    pub integrations: i32,
    /// VirusType enum name (MCV, EBV, HPV, HBV, HHV8), None for "null"
    pub interpretation: Option<String>,
    pub percentage_covered: f64,
    pub mean_coverage: f64,
    pub expected_clonal_coverage: Option<f64>,
    pub reported: bool,
    /// None without the column (older files)
    pub blacklisted: Option<bool>,
    /// VirusLikelihoodType name: HIGH, LOW or UNKNOWN
    pub driver_likelihood: String,
}

/// AnnotatedVirusFile.generateFileName
pub fn annotated_virus_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.virus.annotated.tsv", dir(base_path))
}

/// VirusType.fromVirusName: by the display names (HHV-8 for HHV8)
pub fn virus_type_from_name(name: &str) -> Result<String> {
    Ok(match name {
        "MCV" | "EBV" | "HPV" | "HBV" => name.to_string(),
        "HHV-8" => "HHV8".to_string(),
        _ => bail!("java.lang.IllegalStateException: Cannot resolve virus name: {name}"),
    })
}

/// AnnotatedVirusFile.read
pub fn read_annotated_viruses(path: &str) -> Result<Vec<AnnotatedVirus>> {
    let lines = lines_of(path)?;
    let f = header(&lines, path)?;
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        let expected_clonal_coverage = match f.get("expectedClonalCoverage") {
            Some(&i) => {
                let s = r.at(i)?;
                if s == "null" { None } else { Some(parse_double_value(s)?) }
            }
            None => None,
        };
        let interpretation = {
            let s = r.str("interpretation")?;
            if s == "null" { None } else { Some(virus_type_from_name(s)?) }
        };
        out.push(AnnotatedVirus {
            taxid: r.int("taxid")?,
            name: r.string("name")?,
            qc_status: r.enum_value("qcStatus", &VIRUS_BREAKEND_QC_STATUSES, "VirusBreakendQCStatus")?,
            integrations: r.int("integrations")?,
            interpretation,
            percentage_covered: r.double_or("percentageCovered", 0.0)?,
            mean_coverage: r.double_or("meanCoverage", 0.0)?,
            expected_clonal_coverage,
            reported: r.bool("reported")?,
            blacklisted: match f.get("blacklisted") {
                Some(&i) => Some(crate::java_text::parse_bool(r.at(i)?)),
                None => None,
            },
            driver_likelihood: match f.get("driverLikelihood") {
                Some(&i) => enum_value(r.at(i)?, &["HIGH", "LOW", "UNKNOWN"], "VirusLikelihoodType")?,
                None => "UNKNOWN".to_string(),
            },
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// CUPPA

/// DataType names
pub const CUPPA_DATA_TYPES: [&str; 5] = ["PROB", "FEAT_CONTRIB", "SIG_QUANTILE", "CV_PERFORMANCE", "NONE"];
/// ClassifierGroup names
pub const CUPPA_CLASSIFIER_GROUPS: [&str; 4] = ["COMBINED", "DNA", "RNA", "NONE"];
/// ClassifierName names
pub const CUPPA_CLASSIFIER_NAMES: [&str; 9] =
    ["COMBINED", "DNA_COMBINED", "GEN_POS", "SNV96", "EVENT", "RNA_COMBINED", "GENE_EXP", "ALT_SJ", "NONE"];

#[derive(Clone, Debug, PartialEq)]
pub struct CuppaPredictionEntry {
    pub sample_id: String,
    pub data_type: String,
    pub classifier_group: String,
    pub classifier_name: String,
    pub feature_name: String,
    pub feature_value: f64,
    pub cancer_type: String,
    pub data_value: f64,
    pub rank: i32,
    pub rank_group: i32,
}

/// CuppaPredictions.generateVisDataTsvFilename
pub fn cuppa_vis_data_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.cuppa.vis_data.tsv", dir(base_path))
}

/// CuppaPredictions.generateVisPlotFilename
pub fn cuppa_vis_plot_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.cuppa.vis.png", dir(base_path))
}

/// CuppaPredictions.parseDouble: empty is NaN, "inf" and "-inf" the infinities
fn cuppa_double(s: &str) -> Result<f64> {
    match s {
        "" => Ok(f64::NAN),
        "inf" => Ok(f64::INFINITY),
        "-inf" => Ok(f64::NEG_INFINITY),
        _ => parse_double_value(s),
    }
}

/// CuppaPredictions.parseStringWithEmpty, then toUpperCase
fn cuppa_name(s: &str) -> String {
    if s.is_empty() { "NONE".to_string() } else { s.to_uppercase() }
}

/// CuppaPredictions.fromTsv: the sample-level rows (PROB, FEAT_CONTRIB, SIG_QUANTILE); Err for an empty
/// list, as the CuppaPredictions constructor throws
pub fn read_cuppa_predictions(path: &str) -> Result<Vec<CuppaPredictionEntry>> {
    let lines = lines_of(path)?;
    let f = header(&lines, path)?;
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        // the indices are looked up (unboxed) before the first line
        let data_type = enum_value(&cuppa_name(r.str("data_type")?), &CUPPA_DATA_TYPES, "DataType")?;
        if !matches!(data_type.as_str(), "PROB" | "FEAT_CONTRIB" | "SIG_QUANTILE") {
            continue;
        }
        out.push(CuppaPredictionEntry {
            sample_id: r.string("sample_id")?,
            data_type,
            classifier_group: enum_value(&cuppa_name(r.str("clf_group")?), &CUPPA_CLASSIFIER_GROUPS, "ClassifierGroup")?,
            classifier_name: enum_value(&cuppa_name(r.str("clf_name")?), &CUPPA_CLASSIFIER_NAMES, "ClassifierName")?,
            feature_name: {
                let s = r.str("feat_name")?;
                if s.is_empty() { "NONE".to_string() } else { s.to_string() }
            },
            feature_value: cuppa_double(r.str("feat_value")?)?,
            cancer_type: r.string("cancer_type")?,
            data_value: cuppa_double(r.str("data_value")?)?,
            rank: parse_int_value(r.str("rank")?)?,
            rank_group: parse_int_value(r.str("rank_group")?)?,
        });
    }
    for name in ["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value", "cancer_type", "data_value", "rank", "rank_group"] {
        if !f.contains_key(name) {
            bail!("{path}: java.lang.NullPointerException: no column {name}");
        }
    }
    if out.is_empty() {
        bail!("java.lang.IllegalStateException: `CuppaPredictions` must be instantiated with a non-empty list");
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// PipelineVersionFile

/// PipelineVersionFile.majorDotMinorVersion: Ok(None) for a file of several lines (with a warning) or a
/// version without a dot; Err for an empty or missing file, and where `split("\\.")` leaves no second part
pub fn pipeline_major_dot_minor_version(path: &str) -> Result<Option<String>> {
    let lines = lines_of(path)?;
    if lines.is_empty() {
        bail!("java.io.IOException: Pipeline version file seems empty on {path}");
    }
    if lines.len() > 1 {
        crate::warn!("Too many lines in pipeline version file {path}!");
        return Ok(None);
    }
    let version = &lines[0];
    if !version.contains('.') {
        return Ok(None);
    }
    let parts = split(version, '.');
    if parts.len() < 2 {
        bail!("java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length {}", parts.len());
    }
    Ok(Some(format!("{}.{}", parts[0], parts[1])))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn versions_and_names() {
        let dir = std::env::temp_dir().join(format!("tool_results_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("v.txt");
        std::fs::write(&p, "3.0.0 [oncoanalyser]\n").unwrap();
        assert_eq!(pipeline_major_dot_minor_version(p.to_str().unwrap()).unwrap().as_deref(), Some("3.0"));
        std::fs::write(&p, "5\n").unwrap();
        assert_eq!(pipeline_major_dot_minor_version(p.to_str().unwrap()).unwrap(), None);
        std::fs::write(&p, "5.\n").unwrap();
        assert!(pipeline_major_dot_minor_version(p.to_str().unwrap()).is_err());
        std::fs::write(&p, "").unwrap();
        assert!(pipeline_major_dot_minor_version(p.to_str().unwrap()).is_err());
        assert_eq!(virus_type_from_name("HHV-8").unwrap(), "HHV8");
        assert!(virus_type_from_name("HHV8").is_err());
        let e = snv_signature_etiologies();
        assert_eq!(e.get(&"Sig2".to_string()).map(String::as_str), Some("APOBEC activity"));
        let g = PeachGenotype {
            gene: "G".into(),
            allele: "*1".into(),
            allele_count: 2,
            function: String::new(),
            linked_drugs: String::new(),
            url_prescription_info: String::new(),
            panel_version: String::new(),
            repo_version: String::new(),
        };
        assert_eq!(g.haplotype().unwrap(), "*1_HOM");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
