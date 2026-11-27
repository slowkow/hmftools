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
