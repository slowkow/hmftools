//! hmf-common `rna`: ISOFOX's transcript expression file (TranscriptExpressionFile) and the cohort
//! expression matrix (RnaExpressionMatrix, read through MatrixFile.loadMatrixDataFile).

use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;

use crate::java_text::{fields_index_map, parse_double, read_lines, split_all};
use crate::neo::check_add_dir_separator;

/// RnaCommon.ISF_FILE_ID
pub const ISF_FILE_ID: &str = ".isf.";

/// TranscriptExpression: transcript name, gene name, adjusted TPM
#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptExpression {
    pub transcript_name: String,
    pub gene_name: String,
    pub tpm: f64,
}

/// TranscriptExpressionFile.generateFilename: `<dir>/<sample>.isf.transcript_data.tsv`
pub fn transcript_expression_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{ISF_FILE_ID}transcript_data.tsv", check_add_dir_separator(base_path))
}

/// TranscriptExpressionFile.read: the delimiter follows the extension; TransName, GeneName and AdjTPM
pub fn read_transcript_expression(filename: &str) -> Result<Vec<TranscriptExpression>> {
    let lines = read_lines(filename).with_context(|| format!("java.nio.file.NoSuchFileException: {filename}"))?;
    let header = lines.first().ok_or_else(|| anyhow!("{filename}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let delim = crate::hla::infer_file_delimiter(filename);
    let f = fields_index_map(header, delim);
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("{filename}: missing column {n}"));
    let (trans_i, gene_i, tpm_i) = (col("TransName")?, col("GeneName")?, col("AdjTPM")?);
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = split_all(line, delim);
        let g = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: index {i}"));
        let tpm_s = g(tpm_i)?;
        out.push(TranscriptExpression {
            transcript_name: g(trans_i)?.to_string(),
            gene_name: g(gene_i)?.to_string(),
            tpm: parse_double(tpm_s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{tpm_s}\""))?,
        });
    }
    Ok(out)
}

/// RnaExpressionMatrix.INVALID_EXP
pub const INVALID_EXP: f64 = -1.0;

/// RnaExpressionMatrix: expression by gene id or transcript name and sample, from a CSV matrix with the
/// samples as columns
pub struct RnaExpressionMatrix {
    rows: Vec<Vec<f64>>,
    sample_index: HashMap<String, usize>,
    id_index: HashMap<String, usize>,
}

pub const EXPRESSION_SCOPE_GENE: &str = "GENE";
pub const EXPRESSION_SCOPE_TRANS: &str = "TRANS";

impl RnaExpressionMatrix {
    /// new RnaExpressionMatrix(file, scope). As Java, the matrix is read through a reader that
    /// decompresses `.gz` but the row ids with Files.readAllLines, which does not: Java's matrix is null
    /// for a `.gz` file, an error here.
    pub fn load(filename: &str, scope: &str) -> Result<RnaExpressionMatrix> {
        let transcript_scope = scope == EXPRESSION_SCOPE_TRANS;
        let mut ignore_fields: Vec<&str> = vec!["GeneId", "GeneName"];
        if transcript_scope {
            ignore_fields.push("TransName");
        }
        if filename.ends_with(".gz") {
            anyhow::bail!("failed to read RNA expression data from file({filename}): MalformedInputException");
        }
        let lines = read_lines(filename).with_context(|| format!("failed to read matrix data file({filename})"))?;
        if lines.len() <= 1 {
            anyhow::bail!("empty matrix data");
        }
        // MatrixFile.loadMatrixDataFile
        let header: Vec<&str> = split_all(&lines[0], ',');
        let mut column_names: Vec<String> = header.iter().map(|s| s.to_string()).collect();
        let mut ignore_cols: Vec<usize> = Vec::new();
        for (i, c) in header.iter().enumerate() {
            if ignore_fields.contains(c) {
                ignore_cols.push(i);
                if ignore_cols.len() == ignore_fields.len() {
                    break;
                }
            }
        }
        for f in &ignore_fields {
            if let Some(p) = column_names.iter().position(|c| c == f) {
                column_names.remove(p);
            }
        }
        let col_count = column_names.len();
        let mut rows = vec![vec![0.0; col_count]; lines.len() - 1];
        let mut invalid = 0;
        for (r, line) in lines[1..].iter().enumerate() {
            let items = split_all(line, ',');
            if items.len() != col_count + ignore_cols.len() {
                invalid += 1;
                continue;
            }
            let mut c = 0;
            for (i, item) in items.iter().enumerate() {
                if ignore_cols.contains(&i) {
                    continue;
                }
                let v = if item.is_empty() {
                    0.0
                } else {
                    parse_double(item).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{item}\""))?
                };
                rows[r][c] = v;
                c += 1;
            }
        }
        if invalid > 0 {
            crate::warn!("matrix(rows={} cols={}) has {} invalid rows", rows.len(), col_count, invalid);
        }
        crate::info!("loaded matrix(rows={} cols={}) from file({})", rows.len(), col_count, filename);
        let sample_index: HashMap<String, usize> = column_names.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect();
        let f = fields_index_map(&lines[0], ',');
        let id_col = if transcript_scope { f.get("TransName") } else { f.get("GeneId") }
            .copied()
            .ok_or_else(|| anyhow!("{filename}: missing id column"))?;
        let mut id_index = HashMap::new();
        for (i, line) in lines[1..].iter().enumerate() {
            let items = split_all(line, ',');
            let id = items.get(id_col).ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException"))?;
            id_index.insert(id.to_string(), i);
        }
        Ok(RnaExpressionMatrix { rows, sample_index, id_index })
    }

    /// getExpression: INVALID_EXP for an unknown id or sample
    pub fn expression(&self, gene_trans_id: &str, sample_id: &str) -> f64 {
        match (self.id_index.get(gene_trans_id), self.sample_index.get(sample_id)) {
            (Some(&r), Some(&c)) => self.rows[r][c],
            _ => INVALID_EXP,
        }
    }

    pub fn has_sample_id(&self, sample_id: &str) -> bool {
        self.sample_index.contains_key(sample_id)
    }
    pub fn has_trans_name(&self, name: &str) -> bool {
        self.id_index.contains_key(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_and_transcripts() {
        let dir = std::env::temp_dir().join(format!("rna_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let m = dir.join("m.csv");
        std::fs::write(&m, "GeneId,GeneName,TransName,S1,S2\nG1,N1,T1,1.5,2\nG1,N1,T2,,3\nG2,N2,T3,1\n").unwrap();
        let mat = RnaExpressionMatrix::load(m.to_str().unwrap(), EXPRESSION_SCOPE_TRANS).unwrap();
        assert_eq!(mat.expression("T1", "S2"), 2.0);
        assert_eq!(mat.expression("T2", "S1"), 0.0);
        // a short row is left at zero
        assert_eq!(mat.expression("T3", "S1"), 0.0);
        assert_eq!(mat.expression("T9", "S1"), INVALID_EXP);
        assert!(mat.has_sample_id("S1") && !mat.has_sample_id("GeneId"));
        let t = dir.join("S.isf.transcript_data.tsv");
        std::fs::write(&t, "GeneId\tGeneName\tTransName\tAdjTPM\nG1\tN1\tT1\t4.25\n").unwrap();
        let e = read_transcript_expression(&transcript_expression_filename(dir.to_str().unwrap(), "S")).unwrap();
        assert_eq!(e, vec![TranscriptExpression { transcript_name: "T1".into(), gene_name: "N1".into(), tpm: 4.25 }]);
        std::fs::remove_dir_all(&dir).ok();
    }
}

// ---------------------------------------------------------------------------------------------
// Additions for ORANGE (2026-11-27): RnaStatisticFile (`.isf.summary.tsv`), GeneExpressionFile
// (`.isf.gene_data.tsv`), RnaFusionFile (`.isf.pass_fusions.tsv`) and NovelSpliceJunctionFile
// (`.isf.alt_splice_junc.tsv`), as their readers parse them. Enum values kept as their names.

use crate::linx_files::{enum_value, parse_double_value, parse_int_value, Row};

/// RnaQcFilter names
pub const RNA_QC_FILTERS: [&str; 4] = ["PASS", "FAIL_LOW_COVERAGE", "WARN_DUPLICATE_RATE", "WARN_SPLICED_GENE_COVERAGE"];

/// RnaStatistics (the fields ORANGE reads, and those parsed on the way)
#[derive(Clone, Debug, PartialEq)]
pub struct RnaStatistics {
    pub qc_status: Vec<String>,
    pub total_fragments: i64,
    pub duplicate_fragments: i64,
    pub spliced_fragment_perc: f64,
    pub unspliced_fragment_perc: f64,
    pub alt_fragment_perc: f64,
    pub chimeric_fragment_perc: f64,
    pub spliced_gene_count: i32,
    pub read_length: i32,
}

pub fn rna_statistics_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{ISF_FILE_ID}summary.tsv", check_add_dir_separator(base_path))
}

/// RnaStatisticFile.calcQcStatus
pub fn rna_calc_qc_status(total: i64, duplicates: i64, spliced_genes: i32) -> Vec<String> {
    if total - duplicates < 2_500_000 {
        return vec!["FAIL_LOW_COVERAGE".into()];
    }
    let mut v = Vec::new();
    let dup = if total > 0 { duplicates as f64 / total as f64 } else { 0.0 };
    if dup > 0.9 {
        v.push("WARN_DUPLICATE_RATE".to_string());
    }
    if spliced_genes >= 0 && spliced_genes < 17000 {
        v.push("WARN_SPLICED_GENE_COVERAGE".to_string());
    }
    if v.is_empty() { vec!["PASS".into()] } else { v }
}

/// RnaStatisticFile.fromLines (the delimiter from the header: ',' if it has one; values split without a limit)
pub fn read_rna_statistics(path: &str) -> Result<RnaStatistics> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let header = lines.first().ok_or_else(|| anyhow!("{path}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let delim = if header.contains(',') { ',' } else { '\t' };
    let f = fields_index_map(header, delim);
    let line = lines.get(1).ok_or_else(|| anyhow!("{path}: IndexOutOfBoundsException: Index 1 out of bounds for length 1"))?;
    let r = Row { file: path, fields: &f, values: crate::java_text::split(line, delim) };
    let total = r.long_or("TotalFragments", 0)?;
    let duplicates = r.long_or("DuplicateFragments", 0)?;
    let spliced_genes = r.int_or("SplicedGeneCount", 0)?;
    let qc_status = match f.get("QcStatus") {
        Some(&i) => {
            let s = r.at(i)?;
            s.split(';').map(|x| enum_value(x, &RNA_QC_FILTERS, "RnaQcFilter")).collect::<Result<Vec<_>>>()?
        }
        None => rna_calc_qc_status(total, duplicates, spliced_genes),
    };
    Ok(RnaStatistics {
        qc_status,
        total_fragments: total,
        duplicate_fragments: duplicates,
        spliced_fragment_perc: r.double_or("SplicedFragmentPerc", 0.0)?,
        unspliced_fragment_perc: r.double_or("UnsplicedFragmentPerc", 0.0)?,
        alt_fragment_perc: r.double_or("AltFragmentPerc", 0.0)?,
        chimeric_fragment_perc: r.double_or("ChimericFragmentPerc", 0.0)?,
        spliced_gene_count: r.int_or("SplicedGeneCount", 0)?,
        read_length: {
            let v = r.int_or("ReadLength", 0)?;
            for n in ["FragLength5th", "FragLength50th", "FragLength95th", "MedianGCRatio", "ForwardStrandPercent"] {
                r.double_or(n, 0.0)?;
            }
            v
        },
    })
}

/// GeneExpression (Isofox's gene data)
#[derive(Clone, Debug, PartialEq)]
pub struct GeneExpression {
    pub gene_name: String,
    pub tpm: f64,
    pub spliced_fragments: i32,
    pub unspliced_fragments: i32,
    pub median_tpm_cancer: f64,
    pub percentile_cancer: f64,
    pub median_tpm_cohort: f64,
    pub percentile_cohort: f64,
    pub reported_status: String,
}

pub fn gene_expression_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{ISF_FILE_ID}gene_data.tsv", check_add_dir_separator(base_path))
}

/// GeneExpressionFile.read (Java returns null on an IOException; here an Err)
pub fn read_gene_expressions(path: &str) -> Result<Vec<GeneExpression>> {
    let lines = read_lines(path).with_context(|| format!("failed to load Isofox gene data file({path})"))?;
    let delim = crate::hla::infer_file_delimiter(path);
    let header = lines.first().ok_or_else(|| anyhow!("{path}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let f = fields_index_map(header, delim);
    for n in ["GeneName", "AdjTPM", "SplicedFragments", "UnsplicedFragments"] {
        if !f.contains_key(n) {
            bail_npe(path, n)?;
        }
    }
    let cohort = |r: &Row, n: &str| -> Result<f64> {
        match f.get(n) {
            Some(&i) => parse_double_value(r.at(i)?),
            None => Ok(0.0),
        }
    };
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, delim) };
        let reported_status = match f.get("ReportedStatus") {
            Some(&i) => enum_value(r.at(i)?, &["NONE", "NOT_REPORTED", "REPORTED"], "ReportedStatus")?,
            None => "NONE".to_string(),
        };
        out.push(GeneExpression {
            gene_name: r.string("GeneName")?,
            tpm: r.double("AdjTPM")?,
            spliced_fragments: r.int("SplicedFragments")?,
            unspliced_fragments: r.int("UnsplicedFragments")?,
            median_tpm_cancer: cohort(&r, "MedianTpmCancer")?,
            percentile_cancer: cohort(&r, "PercentileCancer")?,
            median_tpm_cohort: cohort(&r, "MedianTpmCohort")?,
            percentile_cohort: cohort(&r, "PercentileCohort")?,
            reported_status,
        });
    }
    Ok(out)
}

fn bail_npe(path: &str, column: &str) -> Result<()> {
    Err(anyhow!("{path}: java.lang.NullPointerException: no column {column}"))
}

/// RnaFusion (Isofox's passing fusions)
#[derive(Clone, Debug, PartialEq)]
pub struct RnaFusion {
    pub name: String,
    /// KnownFusionType name
    pub known_type: String,
    pub chromosome_up: String,
    pub chromosome_down: String,
    pub position_up: i32,
    pub position_down: i32,
    pub orientation_up: i8,
    pub orientation_down: i8,
    pub junction_type_up: String,
    pub junction_type_down: String,
    pub transcript_up: String,
    pub transcript_down: String,
    pub exon_up: i32,
    pub exon_down: i32,
    pub sv_type: String,
    pub split_fragments: i32,
    pub realigned_frags: i32,
    pub discordant_frags: i32,
    pub depth_up: i32,
    pub depth_down: i32,
    pub max_anchor_length_up: i32,
    pub max_anchor_length_down: i32,
    pub cohort_frequency: i32,
}

pub fn rna_fusion_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{ISF_FILE_ID}pass_fusions.tsv", check_add_dir_separator(base_path))
}

const KNOWN_FUSION_TYPES: [&str; 8] =
    ["NONE", "KNOWN_PAIR", "PROMISCUOUS_5", "PROMISCUOUS_3", "ENHANCER_KNOWN_PAIR", "ENHANCER_PROMISCUOUS", "EXON_DEL_DUP", "PROMISCUOUS_ENHANCER_TARGET"];

/// RnaFusionFile.parseKnownType
fn parse_known_type(s: &str) -> String {
    if KNOWN_FUSION_TYPES.contains(&s) {
        s.to_string()
    } else if s.contains("PROM5") {
        "PROMISCUOUS_5".into()
    } else if s.contains("PROM3") {
        "PROMISCUOUS_3".into()
    } else {
        "NONE".into()
    }
}

/// RnaFusionFile.read, with the pre-v2.0 column names
pub fn read_rna_fusions(path: &str) -> Result<Vec<RnaFusion>> {
    let lines = read_lines(path).with_context(|| format!("failed to load Isofox fusion file({path})"))?;
    let delim = crate::hla::infer_file_delimiter(path);
    let header = lines.first().ok_or_else(|| anyhow!("{path}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let f = fields_index_map(header, delim);
    // getOrDefault(new name, get(old name)): unboxing a null index is a NullPointerException
    let idx = |new: &str, old: &str| -> Result<usize> {
        f.get(new).or_else(|| f.get(old)).copied().ok_or_else(|| anyhow!("{path}: java.lang.NullPointerException: no column {new}"))
    };
    let chr_up = idx("ChromosomeUp", "ChrUp")?;
    let chr_down = idx("ChromosomeDown", "ChrDown")?;
    let pos_up = idx("PositionUp", "PosUp")?;
    let pos_down = idx("PositionDown", "PosDown")?;
    let orient_up = idx("OrientationUp", "OrientUp")?;
    let orient_down = idx("OrientationDown", "OrientDown")?;
    let junc_up = idx("JunctionTypeUp", "JuncTypeUp")?;
    let junc_down = idx("JunctionTypeDown", "JuncTypeDown")?;
    let sv_type = idx("SvType", "SVType")?;
    let depth_up = idx("DepthUp", "CoverageUp")?;
    let depth_down = idx("DepthDown", "CoverageDown")?;
    let cohort = idx("CohortFrequency", "CohortCount")?;
    let gene_up = f.get("GeneNameUp").copied();
    let gene_down = f.get("GeneNameDown").copied();
    let known_type = f.get("KnownType").or_else(|| f.get("KnownFusionType")).copied();
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, delim) };
        let byte = |i: usize| -> Result<i8> {
            let s = r.at(i)?;
            crate::java_text::parse_byte(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
        };
        let transcript_up = r.str_or("TranscriptUp", "")?.to_string();
        let transcript_down = r.str_or("TranscriptDown", "")?.to_string();
        let exon_up = match f.get("ExonUp") {
            Some(&i) => parse_int_value(r.at(i)?)?,
            None => -1,
        };
        let exon_down = match f.get("ExonDown") {
            Some(&i) => parse_int_value(r.at(i)?)?,
            None => -1,
        };
        let name = match (gene_up, gene_down) {
            (Some(u), Some(d)) => format!("{}_{}", r.at(u)?, r.at(d)?),
            _ => r.string("Name")?,
        };
        let known_type = match known_type {
            Some(i) => parse_known_type(r.at(i)?),
            None => "NONE".to_string(),
        };
        out.push(RnaFusion {
            name,
            known_type,
            chromosome_up: r.at(chr_up)?.to_string(),
            chromosome_down: r.at(chr_down)?.to_string(),
            position_up: parse_int_value(r.at(pos_up)?)?,
            position_down: parse_int_value(r.at(pos_down)?)?,
            orientation_up: byte(orient_up)?,
            orientation_down: byte(orient_down)?,
            junction_type_up: r.at(junc_up)?.to_string(),
            junction_type_down: r.at(junc_down)?.to_string(),
            transcript_up,
            transcript_down,
            exon_up,
            exon_down,
            sv_type: enum_value(r.at(sv_type)?, &crate::linx_files::SV_TYPES, "StructuralVariantType")?,
            split_fragments: r.int("SplitFrags")?,
            realigned_frags: r.int("RealignedFrags")?,
            discordant_frags: r.int("DiscordantFrags")?,
            depth_up: parse_int_value(r.at(depth_up)?)?,
            depth_down: parse_int_value(r.at(depth_down)?)?,
            max_anchor_length_up: r.int("MaxAnchorLengthUp")?,
            max_anchor_length_down: r.int("MaxAnchorLengthDown")?,
            cohort_frequency: parse_int_value(r.at(cohort)?)?,
        });
    }
    Ok(out)
}

/// NovelSpliceJunction (Isofox's alternative splice junctions)
#[derive(Clone, Debug, PartialEq)]
pub struct NovelSpliceJunction {
    pub gene_name: String,
    pub chromosome: String,
    pub junction_start: i32,
    pub junction_end: i32,
    pub junction_type: String,
    pub transcript_start: String,
    pub transcript_end: String,
    pub exon_start: i32,
    pub exon_end: i32,
    pub fragment_count: i32,
    pub depth_start: i32,
    pub depth_end: i32,
    pub region_start: String,
    pub region_end: String,
    pub bases_start: String,
    pub bases_end: String,
    pub cohort_frequency: i32,
}

pub fn novel_splice_junction_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{ISF_FILE_ID}alt_splice_junc.tsv", check_add_dir_separator(base_path))
}

pub const ALT_SPLICE_JUNCTION_TYPES: [&str; 11] = [
    "SKIPPED_EXONS",
    "NOVEL_5_PRIME",
    "NOVEL_3_PRIME",
    "NOVEL_EXON",
    "NOVEL_INTRON",
    "MIXED_TRANS",
    "INTRONIC",
    "EXON_INTRON",
    "CIRCULAR",
    "CANONICAL",
    "UNKNOWN",
];
pub const ALT_SPLICE_JUNCTION_CONTEXTS: [&str; 5] = ["SPLICE_JUNC", "EXONIC", "INTRONIC", "MIXED", "UNKNOWN"];

/// NovelSpliceJunctionFile.read
pub fn read_novel_splice_junctions(path: &str) -> Result<Vec<NovelSpliceJunction>> {
    let lines = read_lines(path).with_context(|| format!("failed to load Isofox fusion file({path})"))?;
    let delim = crate::hla::infer_file_delimiter(path);
    let header = lines.first().ok_or_else(|| anyhow!("{path}: java.lang.NullPointerException: empty file"))?;
    let f = fields_index_map(header, delim);
    let mut out = Vec::new();
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, delim) };
        for n in ["GeneName", "Chromosome", "SjStart", "SjEnd", "Type", "FragCount", "DepthStart", "DepthEnd", "RegionStart", "RegionEnd"] {
            r.index(n)?;
        }
        let opt_int = |n: &str, d: i32| -> Result<i32> {
            match f.get(n) {
                Some(&i) => parse_int_value(r.at(i)?),
                None => Ok(d),
            }
        };
        let chromosome = r.string("Chromosome")?;
        let junction_start = r.int("SjStart")?;
        let junction_end = r.int("SjEnd")?;
        out.push(NovelSpliceJunction {
            gene_name: r.string("GeneName")?,
            chromosome,
            junction_start,
            junction_end,
            junction_type: r.enum_value("Type", &ALT_SPLICE_JUNCTION_TYPES, "AltSpliceJunctionType")?,
            transcript_start: r.str_or("TranscriptStart", "")?.to_string(),
            transcript_end: r.str_or("TranscriptEnd", "")?.to_string(),
            exon_start: opt_int("ExonStart", -1)?,
            exon_end: opt_int("ExonEnd", -1)?,
            fragment_count: r.int("FragCount")?,
            depth_start: r.int("DepthStart")?,
            depth_end: r.int("DepthEnd")?,
            region_start: r.enum_value("RegionStart", &ALT_SPLICE_JUNCTION_CONTEXTS, "AltSpliceJunctionContext")?,
            region_end: r.enum_value("RegionEnd", &ALT_SPLICE_JUNCTION_CONTEXTS, "AltSpliceJunctionContext")?,
            bases_start: r.str_or("BaseStart", "")?.to_string(),
            bases_end: r.str_or("BaseEnd", "")?.to_string(),
            cohort_frequency: opt_int("CohortFrequency", 0)?,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod orange_rna_tests {
    use super::*;

    #[test]
    fn qc_status() {
        assert_eq!(rna_calc_qc_status(1_000_000, 0, 20000), vec!["FAIL_LOW_COVERAGE"]);
        assert_eq!(rna_calc_qc_status(3_000_000, 1_000_000, 20000), vec!["FAIL_LOW_COVERAGE"]);
        assert_eq!(rna_calc_qc_status(100_000_000, 95_000_000, 20000), vec!["WARN_DUPLICATE_RATE"]);
        assert_eq!(rna_calc_qc_status(100_000_000, 91_000_000, 0), vec!["WARN_DUPLICATE_RATE", "WARN_SPLICED_GENE_COVERAGE"]);
        assert_eq!(rna_calc_qc_status(100_000_000, 1, 17000), vec!["PASS"]);
    }
}
