//! hmf-common `hla.LilacAllele`: LILAC's allele file `<sample>.lilac.tsv` (read).
//!
//! As LilacAllele.read: an absent `.tsv` falls back to the `.csv` name (checkFileExtensionRename), the
//! delimiter follows the extension (inferFileDelimiter), lines are split without a limit (trailing empty
//! fields dropped) and the Genes column is optional (MHC_CLASS_1 without it).

use anyhow::{anyhow, Context, Result};

use crate::java_text::{fields_index_map, parse_double, parse_int, read_lines, split};

pub const MHC_CLASS_I: &str = "MHC_CLASS_1";
pub const LILAC_FILE_EXTENSION: &str = ".lilac.tsv";

/// LilacAllele
#[derive(Clone, Debug, PartialEq)]
pub struct LilacAllele {
    pub genes: String,
    pub allele: String,
    pub ref_fragments: i32,
    pub ref_unique: i32,
    pub ref_shared: i32,
    pub ref_wild: i32,
    pub tumor_fragments: i32,
    pub tumor_unique: i32,
    pub tumor_shared: i32,
    pub tumor_wild: i32,
    pub rna_fragments: i32,
    pub rna_unique: i32,
    pub rna_shared: i32,
    pub rna_wild: i32,
    pub tumor_copy_number: f64,
    pub somatic_missense: f64,
    pub somatic_nonsense_or_frameshift: f64,
    pub somatic_splice: f64,
    pub somatic_synonymous: f64,
    pub somatic_inframe_indel: f64,
}

/// LilacAllele.generateFilename: `basePath + File.separator + sample + ".lilac.tsv"` (no check for a
/// separator already there)
pub fn lilac_allele_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}{LILAC_FILE_EXTENSION}")
}

/// FileDelimiters.switchFilenameExtension
fn switch_extension(filename: &str, required: &str) -> String {
    if filename.ends_with(required) {
        return filename.to_string();
    }
    let other = if required == ".csv" { ".tsv" } else { ".csv" };
    match filename.strip_suffix(other) {
        Some(stem) => format!("{stem}{required}"),
        None => filename.to_string(),
    }
}

/// FileDelimiters.checkFileExtensionRename
pub fn check_file_extension_rename(filename: &str) -> String {
    if std::path::Path::new(filename).exists() {
        return filename.to_string();
    }
    let switched = if filename.ends_with(".csv") { switch_extension(filename, ".tsv") } else { switch_extension(filename, ".csv") };
    if std::path::Path::new(&switched).exists() { switched } else { filename.to_string() }
}

/// FileDelimiters.inferFileDelimiter
pub fn infer_file_delimiter(filename: &str) -> char {
    if filename.ends_with(".csv") || filename.ends_with(".csv.gz") { ',' } else { '\t' }
}

/// LilacAllele.read
pub fn read_lilac_alleles(path: &str) -> Result<Vec<LilacAllele>> {
    let filename = check_file_extension_rename(path);
    let delim = infer_file_delimiter(&filename);
    let lines = read_lines(&filename).with_context(|| format!("java.nio.file.NoSuchFileException: {filename}"))?;
    let header = lines.first().ok_or_else(|| anyhow!("{filename}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let f = fields_index_map(header, delim);
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("{filename}: missing column {n}"));
    let genes_index = f.get("Genes").copied();
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v = split(line, delim);
        let g = |n: &str| -> Result<&str> {
            let i = col(n)?;
            v.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: index {i} of {}", v.len()))
        };
        let int = |n: &str| -> Result<i32> {
            let s = g(n)?;
            parse_int(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
        };
        let dbl = |n: &str| -> Result<f64> {
            let s = g(n)?;
            parse_double(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
        };
        let allele = g("Allele")?.to_string();
        let genes = match genes_index {
            Some(i) => v.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: index {i}"))?.to_string(),
            None => MHC_CLASS_I.to_string(),
        };
        out.push(LilacAllele {
            genes,
            allele,
            ref_fragments: int("RefTotal")?,
            ref_unique: int("RefUnique")?,
            ref_shared: int("RefShared")?,
            ref_wild: int("RefWild")?,
            tumor_fragments: int("TumorTotal")?,
            tumor_unique: int("TumorUnique")?,
            tumor_shared: int("TumorShared")?,
            tumor_wild: int("TumorWild")?,
            rna_fragments: int("RnaTotal")?,
            rna_unique: int("RnaUnique")?,
            rna_shared: int("RnaShared")?,
            rna_wild: int("RnaWild")?,
            tumor_copy_number: dbl("TumorCopyNumber")?,
            somatic_missense: dbl("SomaticMissense")?,
            somatic_nonsense_or_frameshift: dbl("SomaticNonsenseOrFrameshift")?,
            somatic_splice: dbl("SomaticSplice")?,
            somatic_synonymous: dbl("SomaticSynonymous")?,
            somatic_inframe_indel: dbl("SomaticInframeIndel")?,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_alleles() {
        let dir = std::env::temp_dir().join(format!("lilac_allele_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let header = "Allele\tRefTotal\tRefUnique\tRefShared\tRefWild\tTumorTotal\tTumorUnique\tTumorShared\tTumorWild\tRnaTotal\tRnaUnique\tRnaShared\tRnaWild\tTumorCopyNumber\tSomaticMissense\tSomaticNonsenseOrFrameshift\tSomaticSplice\tSomaticSynonymous\tSomaticInframeIndel";
        let line = "A*01:01\t10\t8\t2\t0\t12\t9\t3\t0\t0\t0\t0\t0\t1.5\t1.0\t0.0\t0.0\t0.0\t0.0";
        // written as .csv with tabs: the reader falls back to the .csv name and splits at commas
        let path = dir.join("S.lilac.tsv");
        std::fs::write(&path, format!("{header}\n{line}\n")).unwrap();
        let a = read_lilac_alleles(&lilac_allele_filename(dir.to_str().unwrap(), "S")).unwrap();
        assert_eq!(a.len(), 1);
        assert_eq!(a[0].genes, MHC_CLASS_I);
        assert_eq!(a[0].allele, "A*01:01");
        assert_eq!(a[0].tumor_copy_number, 1.5);
        assert_eq!(a[0].somatic_missense, 1.0);
        assert_eq!(switch_extension("x.lilac.tsv", ".csv"), "x.lilac.csv");
        std::fs::remove_dir_all(&dir).ok();
    }
}

// ---------------------------------------------------------------------------------------------
// Addition for ORANGE (2026-11-27): LilacQcData, LILAC's `<sample>.lilac.qc.tsv` (read)

/// LilacQcData (the fields the reader parses)
#[derive(Clone, Debug, PartialEq)]
pub struct LilacQcData {
    pub genes: String,
    pub status: String,
    pub total_fragments: i32,
    pub fitted_fragments: i32,
    pub discarded_indels: i32,
    pub discarded_alignment_fragments: i32,
    pub hla_y_allele: String,
}

/// LilacQcData.generateFilename: checkAddDirSeparator(basePath) + sample + ".lilac.qc.tsv"
pub fn lilac_qc_filename(base_path: &str, sample: &str) -> String {
    let base = if base_path.ends_with('/') { base_path.to_string() } else { format!("{base_path}/") };
    format!("{base}{sample}.lilac.qc.tsv")
}

/// LilacQcData.read: as LilacAllele.read, the `.csv` fallback, the delimiter by extension, lines split
/// without a limit and Genes optional
pub fn read_lilac_qc(path: &str) -> Result<Vec<LilacQcData>> {
    let filename = check_file_extension_rename(path);
    let delim = infer_file_delimiter(&filename);
    let lines = read_lines(&filename).with_context(|| format!("java.nio.file.NoSuchFileException: {filename}"))?;
    let header = lines.first().ok_or_else(|| anyhow!("{filename}: IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let f = fields_index_map(header, delim);
    let genes_index = f.get("Genes").copied();
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v = split(line, delim);
        let at = |i: usize| -> Result<&str> {
            v.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: index {i} of {}", v.len()))
        };
        let g = |n: &str| -> Result<&str> { at(f.get(n).copied().ok_or_else(|| anyhow!("{filename}: missing column {n}"))?) };
        let int = |n: &str| -> Result<i32> {
            let s = g(n)?;
            parse_int(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
        };
        let genes = match genes_index {
            Some(i) => at(i)?.to_string(),
            None => MHC_CLASS_I.to_string(),
        };
        out.push(LilacQcData {
            genes,
            status: g("Status")?.to_string(),
            total_fragments: int("TotalFragments")?,
            fitted_fragments: int("FittedFragments")?,
            discarded_alignment_fragments: int("DiscardedAlignmentFragments")?,
            discarded_indels: int("DiscardedIndels")?,
            hla_y_allele: g("HlaYAllele")?.to_string(),
        });
    }
    Ok(out)
}
