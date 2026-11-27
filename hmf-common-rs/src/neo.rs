//! hmf-common `neo`: the neoepitope files passed between LINX, ISOFOX and NEO.
//!
//! - `NeoEpitopeType`;
//! - `NeoEpitopeFile`: NEO's `<sample>.neo.neo_data.tsv` (header, lines, reading);
//! - `NeoEpitopeFusion`: LINX's `<sample>.linx.neoepitope.tsv` (reading);
//! - `RnaNeoEpitope`: ISOFOX's `<sample>.isf.neoepitope.tsv` (reading);
//! - `pointMutationInfo`, `fusionInfo` and `extractTranscriptNames`.
//!
//! The readers split lines as the Java readers do: NeoEpitopeFusion without a limit (trailing empty
//! fields dropped, so a line ending in empty fields is too short, an ArrayIndexOutOfBoundsException in
//! Java and an error here), the others with -1.

use anyhow::{anyhow, Context, Result};

use crate::java_text::{fields_index_map, parse_byte, parse_double, parse_int, read_lines, split, split_all};

pub const NEO_FILE_ID: &str = "neo";
pub const ITEM_DELIM: char = ';';

/// NeoEpitopeType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NeoEpitopeType {
    Missense,
    InframeInsertion,
    InframeDeletion,
    Frameshift,
    StopLost,
    InframeFusion,
    OutOfFrameFusion,
}

impl NeoEpitopeType {
    pub fn name(self) -> &'static str {
        match self {
            NeoEpitopeType::Missense => "MISSENSE",
            NeoEpitopeType::InframeInsertion => "INFRAME_INSERTION",
            NeoEpitopeType::InframeDeletion => "INFRAME_DELETION",
            NeoEpitopeType::Frameshift => "FRAMESHIFT",
            NeoEpitopeType::StopLost => "STOP_LOST",
            NeoEpitopeType::InframeFusion => "INFRAME_FUSION",
            NeoEpitopeType::OutOfFrameFusion => "OUT_OF_FRAME_FUSION",
        }
    }
    /// valueOf: None where Java throws IllegalArgumentException
    pub fn from_name(s: &str) -> Option<NeoEpitopeType> {
        Some(match s {
            "MISSENSE" => NeoEpitopeType::Missense,
            "INFRAME_INSERTION" => NeoEpitopeType::InframeInsertion,
            "INFRAME_DELETION" => NeoEpitopeType::InframeDeletion,
            "FRAMESHIFT" => NeoEpitopeType::Frameshift,
            "STOP_LOST" => NeoEpitopeType::StopLost,
            "INFRAME_FUSION" => NeoEpitopeType::InframeFusion,
            "OUT_OF_FRAME_FUSION" => NeoEpitopeType::OutOfFrameFusion,
            _ => return None,
        })
    }
    pub fn is_fusion(self) -> bool {
        self == NeoEpitopeType::InframeFusion || self == NeoEpitopeType::OutOfFrameFusion
    }
    pub fn is_point_mutation(self) -> bool {
        !self.is_fusion()
    }
}

pub fn parse_neo_epitope_type(s: &str) -> Result<NeoEpitopeType> {
    NeoEpitopeType::from_name(s).ok_or_else(|| anyhow!("No enum constant com.hartwig.hmftools.common.neo.NeoEpitopeType.{s}"))
}

/// FileWriterUtils.checkAddDirSeparator: null and "" are returned as they are
pub fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// NeoEpitopeFile.pointMutationInfo
pub fn point_mutation_info(chromosome: &str, position: i32, ref_bases: &str, alt: &str) -> String {
    format!("{chromosome}:{position}:{ref_bases}:{alt}")
}

/// NeoEpitopeFile.fusionInfo
pub fn fusion_info(chromosomes: &[String; 2], positions: &[i32; 2], orientations: &[i8; 2]) -> String {
    format!(
        "{}:{}:{};{}:{}:{}",
        chromosomes[0], positions[0], orientations[0], chromosomes[1], positions[1], orientations[1]
    )
}

/// NeoEpitopeFile.extractTranscriptNames: String.split(";") without a limit (an empty string gives one
/// empty name)
pub fn extract_transcript_names(trans_up: &str, trans_down: &str) -> (Vec<String>, Vec<String>) {
    (
        split(trans_up, ITEM_DELIM).into_iter().map(String::from).collect(),
        split(trans_down, ITEM_DELIM).into_iter().map(String::from).collect(),
    )
}

// ---------------------------------------------------------------------------------------------
// NeoEpitopeFile

/// NeoEpitopeFile: one line of `<sample>.neo.neo_data.tsv`
#[derive(Clone, Debug, PartialEq)]
pub struct NeoEpitopeFile {
    pub id: i32,
    pub variant_type: NeoEpitopeType,
    pub variant_info: String,
    pub variant_copy_number: f64,
    pub copy_number: f64,
    pub subclonal_likelihood: f64,
    pub gene_ids: [String; 2],
    pub gene_names: [String; 2],
    pub chromosomes: [String; 2],
    pub orientations: [i8; 2],
    pub nmd_bases: [i32; 2],
    pub coding_bases_length: [i32; 2],
    pub upstream_aa: String,
    pub downstream_aa: String,
    pub novel_aa: String,
    pub transcripts: [String; 2],
    pub wildtype_aa: String,
    pub fused_intron_length: i32,
    /// skipped donors (up) and acceptors (down)
    pub skipped_acceptors_donors: [i32; 2],
    pub coding_base_positions: [[i32; 2]; 2],
    pub coding_bases: [String; 2],
    pub coding_base_cigars: [String; 2],
}

pub const NEO_DATA_FILE_EXTENSION: &str = ".neo.neo_data.tsv";

/// NeoEpitopeFile.generateFilename
pub fn neo_epitope_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{NEO_DATA_FILE_EXTENSION}", check_add_dir_separator(base_path))
}

pub const FLD_NE_ID: &str = "NeId";
pub const FLD_NE_VAR_TYPE: &str = "VariantType";
pub const FLD_NE_VAR_INFO: &str = "VariantInfo";
pub const FLD_NE_GENE_ID_UP: &str = "GeneIdUp";
pub const FLD_NE_GENE_ID_DOWN: &str = "GeneIdDown";
pub const FLD_NE_GENE_NAME_UP: &str = "GeneNameUp";
pub const FLD_NE_GENE_NAME_DOWN: &str = "GeneNameDown";
pub const FLD_NE_AA_UP: &str = "UpstreamAA";
pub const FLD_NE_AA_DOWN: &str = "DownstreamAA";
pub const FLD_NE_AA_NOVEL: &str = "NovelAA";
pub const FLD_NE_TRANS_UP: &str = "UpTranscripts";
pub const FLD_NE_TRANS_DOWN: &str = "DownTranscripts";
pub const FLD_NE_NMD_MIN: &str = "NmdMin";
pub const FLD_NE_NMD_MAX: &str = "NmdMax";
pub const FLD_NE_CB_LEN_MIN: &str = "CodingBasesLengthMin";
pub const FLD_NE_CB_LEN_MAX: &str = "CodingBasesLengthMax";
pub const FLD_NE_FUSED_LEN: &str = "FusedIntronLength";
pub const FLD_NE_SKIP_DONORS: &str = "SkippedDonors";
pub const FLD_NE_SKIP_ACCEPTORS: &str = "SkippedAcceptors";
pub const FLD_NE_VAR_CN: &str = "VariantCopyNumber";
pub const FLD_NE_CN: &str = "CopyNumber";
pub const FLD_NE_SC_LIKELIHOOD: &str = "SubclonalLikelihood";

/// NeoEpitopeFile.header
pub fn neo_epitope_header() -> String {
    [
        FLD_NE_ID,
        FLD_NE_VAR_TYPE,
        FLD_NE_VAR_INFO,
        FLD_NE_VAR_CN,
        FLD_NE_CN,
        FLD_NE_SC_LIKELIHOOD,
        FLD_NE_GENE_ID_UP,
        FLD_NE_GENE_ID_DOWN,
        FLD_NE_GENE_NAME_UP,
        FLD_NE_GENE_NAME_DOWN,
        "ChrUp",
        "ChrDown",
        "OrientUp",
        "OrientDown",
        FLD_NE_AA_UP,
        FLD_NE_AA_DOWN,
        FLD_NE_AA_NOVEL,
        FLD_NE_NMD_MIN,
        FLD_NE_NMD_MAX,
        FLD_NE_CB_LEN_MIN,
        FLD_NE_CB_LEN_MAX,
        FLD_NE_FUSED_LEN,
        FLD_NE_SKIP_DONORS,
        FLD_NE_SKIP_ACCEPTORS,
        FLD_NE_TRANS_UP,
        FLD_NE_TRANS_DOWN,
        "WildtypeAA",
        "CodingBaseUpPosStart",
        "CodingBaseUpPosEnd",
        "CodingBasesUp",
        "CodingBaseCigarUp",
        "CodingBaseDownPosStart",
        "CodingBaseDownPosEnd",
        "CodingBasesDown",
        "CodingBaseCigarDown",
    ]
    .join("\t")
}

/// NeoEpitopeFile.toString
pub fn neo_epitope_line(neo: &NeoEpitopeFile) -> String {
    use crate::fmt::java_fixed;
    let mut f: Vec<String> = Vec::with_capacity(35);
    f.push(neo.id.to_string());
    f.push(neo.variant_type.name().to_string());
    f.push(neo.variant_info.clone());
    f.push(java_fixed(neo.variant_copy_number, 4));
    f.push(java_fixed(neo.copy_number, 4));
    f.push(java_fixed(neo.subclonal_likelihood, 4));
    f.push(neo.gene_ids[0].clone());
    f.push(neo.gene_ids[1].clone());
    f.push(neo.gene_names[0].clone());
    f.push(neo.gene_names[1].clone());
    f.push(neo.chromosomes[0].clone());
    f.push(neo.chromosomes[1].clone());
    f.push(neo.orientations[0].to_string());
    f.push(neo.orientations[1].to_string());
    f.push(neo.upstream_aa.clone());
    f.push(neo.downstream_aa.clone());
    f.push(neo.novel_aa.clone());
    f.push(neo.nmd_bases[0].to_string());
    f.push(neo.nmd_bases[1].to_string());
    f.push(neo.coding_bases_length[0].to_string());
    f.push(neo.coding_bases_length[1].to_string());
    f.push(neo.fused_intron_length.to_string());
    f.push(neo.skipped_acceptors_donors[0].to_string());
    f.push(neo.skipped_acceptors_donors[1].to_string());
    f.push(neo.transcripts[0].clone());
    f.push(neo.transcripts[1].clone());
    f.push(neo.wildtype_aa.clone());
    f.push(neo.coding_base_positions[0][0].to_string());
    f.push(neo.coding_base_positions[0][1].to_string());
    f.push(neo.coding_bases[0].clone());
    f.push(neo.coding_base_cigars[0].clone());
    f.push(neo.coding_base_positions[1][0].to_string());
    f.push(neo.coding_base_positions[1][1].to_string());
    f.push(neo.coding_bases[1].clone());
    f.push(neo.coding_base_cigars[1].clone());
    f.join("\t")
}

fn col(f: &std::collections::HashMap<String, usize>, name: &str, file: &str) -> Result<usize> {
    // Java unboxes the null index of a missing column: NullPointerException
    f.get(name).copied().ok_or_else(|| anyhow!("{file}: missing column {name}"))
}

fn field<'a>(v: &[&'a str], i: usize, file: &str) -> Result<&'a str> {
    v.get(i).copied().ok_or_else(|| anyhow!("{file}: ArrayIndexOutOfBoundsException: index {i} of {}", v.len()))
}

fn int_of(s: &str) -> Result<i32> {
    parse_int(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

fn dbl_of(s: &str) -> Result<f64> {
    parse_double(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

fn byte_of(s: &str) -> Result<i8> {
    parse_byte(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

/// NeoEpitopeFile.read
pub fn read_neo_epitope_file(path: &str) -> Result<Vec<NeoEpitopeFile>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let mut out = Vec::new();
    if lines.is_empty() {
        return Ok(out);
    }
    let f = fields_index_map(&lines[0], '\t');
    let c = |n: &str| col(&f, n, path);
    let idx = [
        c(FLD_NE_ID)?,
        c(FLD_NE_VAR_TYPE)?,
        c(FLD_NE_VAR_INFO)?,
        c(FLD_NE_VAR_CN)?,
        c(FLD_NE_CN)?,
        c(FLD_NE_SC_LIKELIHOOD)?,
        c(FLD_NE_GENE_ID_UP)?,
        c(FLD_NE_GENE_ID_DOWN)?,
        c(FLD_NE_GENE_NAME_UP)?,
        c(FLD_NE_GENE_NAME_DOWN)?,
        c("ChrUp")?,
        c("ChrDown")?,
        c("OrientUp")?,
        c("OrientDown")?,
        c(FLD_NE_AA_UP)?,
        c(FLD_NE_AA_DOWN)?,
        c(FLD_NE_AA_NOVEL)?,
        c(FLD_NE_NMD_MIN)?,
        c(FLD_NE_NMD_MAX)?,
        c(FLD_NE_CB_LEN_MIN)?,
        c(FLD_NE_CB_LEN_MAX)?,
        c(FLD_NE_FUSED_LEN)?,
        c(FLD_NE_SKIP_DONORS)?,
        c(FLD_NE_SKIP_ACCEPTORS)?,
        c(FLD_NE_TRANS_UP)?,
        c(FLD_NE_TRANS_DOWN)?,
        c("WildtypeAA")?,
        c("CodingBaseUpPosStart")?,
        c("CodingBaseUpPosEnd")?,
        c("CodingBasesUp")?,
        c("CodingBaseCigarUp")?,
        c("CodingBaseDownPosStart")?,
        c("CodingBaseDownPosEnd")?,
        c("CodingBasesDown")?,
        c("CodingBaseCigarDown")?,
    ];
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        let g = |k: usize| field(&v, idx[k], path);
        out.push(NeoEpitopeFile {
            id: int_of(g(0)?)?,
            variant_type: parse_neo_epitope_type(g(1)?)?,
            variant_info: g(2)?.to_string(),
            variant_copy_number: dbl_of(g(3)?)?,
            copy_number: dbl_of(g(4)?)?,
            subclonal_likelihood: dbl_of(g(5)?)?,
            gene_ids: [g(6)?.to_string(), g(7)?.to_string()],
            gene_names: [g(8)?.to_string(), g(9)?.to_string()],
            chromosomes: [g(10)?.to_string(), g(11)?.to_string()],
            orientations: [byte_of(g(12)?)?, byte_of(g(13)?)?],
            upstream_aa: g(14)?.to_string(),
            downstream_aa: g(15)?.to_string(),
            novel_aa: g(16)?.to_string(),
            nmd_bases: [int_of(g(17)?)?, int_of(g(18)?)?],
            coding_bases_length: [int_of(g(19)?)?, int_of(g(20)?)?],
            fused_intron_length: int_of(g(21)?)?,
            skipped_acceptors_donors: [int_of(g(22)?)?, int_of(g(23)?)?],
            transcripts: [g(24)?.to_string(), g(25)?.to_string()],
            wildtype_aa: g(26)?.to_string(),
            coding_base_positions: [[int_of(g(27)?)?, int_of(g(28)?)?], [int_of(g(31)?)?, int_of(g(32)?)?]],
            coding_bases: [g(29)?.to_string(), g(33)?.to_string()],
            coding_base_cigars: [g(30)?.to_string(), g(34)?.to_string()],
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// NeoEpitopeFusion (LINX)

/// NeoEpitopeFusion: one line of LINX's `<sample>.linx.neoepitope.tsv`
#[derive(Clone, Debug, PartialEq)]
pub struct NeoEpitopeFusion {
    pub gene_ids: [String; 2],
    pub gene_names: [String; 2],
    pub chromosomes: [String; 2],
    pub positions: [i32; 2],
    pub orientations: [i8; 2],
    pub sv_ids: [i32; 2],
    pub junction_copy_number: f64,
    pub copy_number: f64,
    pub insert_sequence: String,
    pub chain_length: i32,
    pub transcripts: [String; 2],
}

/// NeoEpitopeFusion.generateFilename
pub fn neo_epitope_fusion_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.linx.neoepitope.tsv", check_add_dir_separator(base_path))
}

/// NeoEpitopeFusion.header
pub fn neo_epitope_fusion_header() -> String {
    [
        "geneIdUp",
        "geneNameUp",
        "chromosomeUp",
        "positionUp",
        "orientationUp",
        "svIdUp",
        "geneIdDown",
        "geneNameDown",
        "chromosomeDown",
        "positionDown",
        "orientationDown",
        "svIdDown",
        "junctionCopyNumber",
        "copyNumber",
        "insertSeq",
        "chainLength",
        "transcriptsUp",
        "transcriptsDown",
    ]
    .join("\t")
}

/// NeoEpitopeFusion.read
pub fn read_neo_epitope_fusions(path: &str) -> Result<Vec<NeoEpitopeFusion>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let mut out = Vec::new();
    if lines.is_empty() {
        return Ok(out);
    }
    let f = fields_index_map(&lines[0], '\t');
    let names = [
        "geneIdUp",
        "geneIdDown",
        "geneNameUp",
        "geneNameDown",
        "chromosomeUp",
        "chromosomeDown",
        "positionUp",
        "positionDown",
        "orientationUp",
        "orientationDown",
        "svIdUp",
        "svIdDown",
        "junctionCopyNumber",
        "copyNumber",
        "insertSeq",
        "chainLength",
        "transcriptsUp",
        "transcriptsDown",
    ];
    let mut idx = [0usize; 18];
    for (i, n) in names.iter().enumerate() {
        idx[i] = col(&f, n, path)?;
    }
    for line in &lines[1..] {
        // line.split(TSV_DELIM): trailing empty fields are dropped
        let v = split(line, '\t');
        let g = |k: usize| field(&v, idx[k], path);
        out.push(NeoEpitopeFusion {
            gene_ids: [g(0)?.to_string(), g(1)?.to_string()],
            gene_names: [g(2)?.to_string(), g(3)?.to_string()],
            chromosomes: [g(4)?.to_string(), g(5)?.to_string()],
            positions: [int_of(g(6)?)?, int_of(g(7)?)?],
            orientations: [byte_of(g(8)?)?, byte_of(g(9)?)?],
            sv_ids: [int_of(g(10)?)?, int_of(g(11)?)?],
            junction_copy_number: dbl_of(g(12)?)?,
            copy_number: dbl_of(g(13)?)?,
            insert_sequence: g(14)?.to_string(),
            chain_length: int_of(g(15)?)?,
            transcripts: [g(16)?.to_string(), g(17)?.to_string()],
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// RnaNeoEpitope (ISOFOX)

/// RnaNeoEpitope: one line of ISOFOX's `<sample>.isf.neoepitope.tsv`
#[derive(Clone, Debug, PartialEq)]
pub struct RnaNeoEpitope {
    pub id: i32,
    pub variant_type: NeoEpitopeType,
    pub variant_info: String,
    pub fragment_count: i32,
    pub base_depth: [i32; 2],
}

/// RnaNeoEpitope.generateFilename: `<dir>/<sample>.isf.neoepitope.tsv`
pub fn rna_neo_epitope_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.isf.neoepitope.tsv", check_add_dir_separator(base_path))
}

/// RnaNeoEpitope.header
pub fn rna_neo_epitope_header() -> String {
    [FLD_NE_ID, FLD_NE_VAR_TYPE, FLD_NE_VAR_INFO, "FragmentCount", "BaseDepthUp", "BaseDepthDown"].join("\t")
}

/// RnaNeoEpitope.read
pub fn read_rna_neo_epitopes(path: &str) -> Result<Vec<RnaNeoEpitope>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let mut out = Vec::new();
    if lines.is_empty() {
        return Ok(out);
    }
    let f = fields_index_map(&lines[0], '\t');
    let idx = [
        col(&f, FLD_NE_ID, path)?,
        col(&f, FLD_NE_VAR_TYPE, path)?,
        col(&f, FLD_NE_VAR_INFO, path)?,
        col(&f, "FragmentCount", path)?,
        col(&f, "BaseDepthUp", path)?,
        col(&f, "BaseDepthDown", path)?,
    ];
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        let g = |k: usize| field(&v, idx[k], path);
        out.push(RnaNeoEpitope {
            id: int_of(g(0)?)?,
            variant_type: parse_neo_epitope_type(g(1)?)?,
            variant_info: g(2)?.to_string(),
            fragment_count: int_of(g(3)?)?,
            base_depth: [int_of(g(4)?)?, int_of(g(5)?)?],
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_and_infos() {
        for t in [
            NeoEpitopeType::Missense,
            NeoEpitopeType::InframeInsertion,
            NeoEpitopeType::InframeDeletion,
            NeoEpitopeType::Frameshift,
            NeoEpitopeType::StopLost,
            NeoEpitopeType::InframeFusion,
            NeoEpitopeType::OutOfFrameFusion,
        ] {
            assert_eq!(NeoEpitopeType::from_name(t.name()), Some(t));
        }
        assert!(NeoEpitopeType::InframeFusion.is_fusion());
        assert!(NeoEpitopeType::StopLost.is_point_mutation());
        assert_eq!(point_mutation_info("1", 100, "A", "T"), "1:100:A:T");
        assert_eq!(fusion_info(&["1".into(), "2".into()], &[100, 200], &[1, -1]), "1:100:1;2:200:-1");
        let (up, down) = extract_transcript_names("T1;T2", "");
        assert_eq!(up, vec!["T1", "T2"]);
        assert_eq!(down, vec![""]);
        assert_eq!(check_add_dir_separator(""), "");
        assert_eq!(check_add_dir_separator("a"), "a/");
    }

    #[test]
    fn neo_file_round_trip() {
        let neo = NeoEpitopeFile {
            id: 3,
            variant_type: NeoEpitopeType::Frameshift,
            variant_info: point_mutation_info("1", 100, "A", "AT"),
            variant_copy_number: 1.23456,
            copy_number: 2.0,
            subclonal_likelihood: 0.0,
            gene_ids: ["G1".into(), "G1".into()],
            gene_names: ["GENE".into(), "GENE".into()],
            chromosomes: ["1".into(), "1".into()],
            orientations: [1, -1],
            nmd_bases: [-1, -1],
            coding_bases_length: [10, 10],
            upstream_aa: "ABC".into(),
            downstream_aa: "".into(),
            novel_aa: "DE_".into(),
            transcripts: ["T1;T2".into(), "T1;T2".into()],
            wildtype_aa: "".into(),
            fused_intron_length: 0,
            skipped_acceptors_donors: [0, 0],
            coding_base_positions: [[90, 120], [0, 0]],
            coding_bases: ["ACGT".into(), "".into()],
            coding_base_cigars: ["4M".into(), "".into()],
        };
        let dir = std::env::temp_dir().join(format!("neo_file_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("S.neo.neo_data.tsv");
        std::fs::write(&path, format!("{}\n{}\n", neo_epitope_header(), neo_epitope_line(&neo))).unwrap();
        let read = read_neo_epitope_file(path.to_str().unwrap()).unwrap();
        assert_eq!(read.len(), 1);
        assert_eq!(read[0].variant_copy_number, 1.2346);
        assert_eq!(read[0].coding_bases, neo.coding_bases);
        assert_eq!(read[0].transcripts, neo.transcripts);
        std::fs::remove_dir_all(&dir).ok();
    }
}
