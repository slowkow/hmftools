//! hmf-common `linx` files as their readers parse them (ORANGE): LinxSvAnnotation (`.linx.svs.tsv`),
//! LinxBreakend (`.linx.breakend.tsv`), LinxCluster (`.linx.clusters.tsv`), LinxFusion (`.linx.fusion.tsv`),
//! LinxDriver (`.linx.drivers.tsv`) and LinxGermlineDisruption (`.linx.germline.disruption.tsv`), with the
//! file names and the breakend coordinate helpers.
//!
//! As Java: the svs, breakend, cluster and fusion lines are split without a limit (`String.split("\t")`, so
//! trailing empty fields are dropped and a column after them is out of range), the driver and disruption
//! lines with -1; enum-valued fields are checked against the enum's names (valueOf) and kept as strings.
//! Every error where Java throws is an `Err`.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

use crate::java_text::{fields_index_map, parse_bool, parse_byte, parse_double, parse_int, read_lines, split, split_all};

/// StructuralVariantType names
pub const SV_TYPES: [&str; 7] = ["BND", "DEL", "DUP", "INF", "INS", "INV", "SGL"];
/// purple.ReportedStatus names
pub const REPORTED_STATUSES: [&str; 3] = ["NONE", "NOT_REPORTED", "REPORTED"];
/// gene.TranscriptRegionType names
pub const TRANSCRIPT_REGION_TYPES: [&str; 6] = ["UNKNOWN", "UPSTREAM", "EXONIC", "INTRONIC", "IG", "DOWNSTREAM"];
/// gene.TranscriptCodingType names
pub const TRANSCRIPT_CODING_TYPES: [&str; 6] = ["UNKNOWN", "CODING", "UTR_5P", "UTR_3P", "NON_CODING", "ENHANCER"];
/// linx.FusionPhasedType names
pub const FUSION_PHASED_TYPES: [&str; 3] = ["INFRAME", "SKIPPED_EXONS", "OUT_OF_FRAME"];
/// linx.FusionLikelihoodType names
pub const FUSION_LIKELIHOOD_TYPES: [&str; 3] = ["HIGH", "LOW", "NA"];
/// linx.FusionReportableReason names
pub const FUSION_REPORTABLE_REASONS: [&str; 15] = [
    "OK",
    "NOT_KNOWN",
    "UNPHASED_NOT_KNOWN",
    "UNPHASED_5P_UTR",
    "UNPHASED_SHORT",
    "SGL_NOT_KNOWN",
    "PRE_GENE_DISTANCE",
    "NMD",
    "NEG_SPLICE_ACC_DISTANCE",
    "EXON_SKIPPING",
    "CHAIN_TERMINATED",
    "NON_DISRUPTIVE_CHAIN",
    "INVALID_TRAVERSAL",
    "CHAIN_LINKS",
    "PROTEIN_DOMAINS",
];
/// linx.DriverEventType names
pub const DRIVER_EVENT_TYPES: [&str; 12] = [
    "UNCLEAR",
    "GAIN",
    "GAIN_ARM",
    "GAIN_CHR",
    "LOH",
    "LOH_ARM",
    "LOH_CHR",
    "LOH_SV_CENTRO",
    "LOH_SV_TELO",
    "DEL",
    "HOM_DEL_DISRUPTION",
    "HOM_DUP_DISRUPTION",
];

/// LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM
pub const BREAKEND_ORIENTATION_UPSTREAM: &str = "Upstream";
/// LinxCommonTypes.SV_VIS_CLUSTER_PREFIX
pub const SV_VIS_CLUSTER_PREFIX: &str = ".cluster";

/// LinxCommonTypes.formVisPlotClusterPrefix: ".cluster-<id>"
pub fn vis_plot_cluster_prefix(cluster_id: &str) -> String {
    format!("{SV_VIS_CLUSTER_PREFIX}-{cluster_id}")
}

/// Enum.valueOf: the name itself when it is one of `names`
pub fn enum_value(value: &str, names: &[&str], enum_name: &str) -> Result<String> {
    if names.contains(&value) {
        Ok(value.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant {enum_name}.{value}")
    }
}

/// A data line of a TSV with Java's field access: a missing column is a NullPointerException (unboxing the
/// null index), a short line an ArrayIndexOutOfBoundsException, a bad number a NumberFormatException.
pub struct Row<'a> {
    pub file: &'a str,
    pub fields: &'a HashMap<String, usize>,
    pub values: Vec<&'a str>,
}

impl<'a> Row<'a> {
    pub fn index(&self, name: &str) -> Result<usize> {
        self.fields.get(name).copied().ok_or_else(|| anyhow!("{}: java.lang.NullPointerException: no column {name}", self.file))
    }
    pub fn at(&self, i: usize) -> Result<&'a str> {
        self.values.get(i).copied().ok_or_else(|| {
            anyhow!("{}: java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", self.file, self.values.len())
        })
    }
    pub fn str(&self, name: &str) -> Result<&'a str> {
        self.at(self.index(name)?)
    }
    pub fn string(&self, name: &str) -> Result<String> {
        Ok(self.str(name)?.to_string())
    }
    /// an optional column: `default` without it
    pub fn str_or(&self, name: &str, default: &'a str) -> Result<&'a str> {
        match self.fields.get(name) {
            Some(&i) => self.at(i),
            None => Ok(default),
        }
    }
    pub fn int(&self, name: &str) -> Result<i32> {
        parse_int_value(self.str(name)?)
    }
    pub fn long(&self, name: &str) -> Result<i64> {
        let s = self.str(name)?;
        crate::java_text::parse_long(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
    }
    pub fn double(&self, name: &str) -> Result<f64> {
        parse_double_value(self.str(name)?)
    }
    pub fn bool(&self, name: &str) -> Result<bool> {
        Ok(parse_bool(self.str(name)?))
    }
    pub fn byte(&self, name: &str) -> Result<i8> {
        let s = self.str(name)?;
        parse_byte(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: Value out of range. Value:\"{s}\""))
    }
    /// FileReaderUtils.getIntValue(map, field, default, values): the default without the column
    pub fn int_or(&self, name: &str, default: i32) -> Result<i32> {
        match self.fields.get(name) {
            Some(&i) => parse_int_value(self.at(i)?),
            None => Ok(default),
        }
    }
    pub fn long_or(&self, name: &str, default: i64) -> Result<i64> {
        match self.fields.get(name) {
            Some(&i) => {
                let s = self.at(i)?;
                crate::java_text::parse_long(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
            }
            None => Ok(default),
        }
    }
    pub fn double_or(&self, name: &str, default: f64) -> Result<f64> {
        match self.fields.get(name) {
            Some(&i) => parse_double_value(self.at(i)?),
            None => Ok(default),
        }
    }
    pub fn enum_value(&self, name: &str, names: &[&str], enum_name: &str) -> Result<String> {
        enum_value(self.str(name)?, names, enum_name)
    }
}

pub fn parse_int_value(s: &str) -> Result<i32> {
    parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

pub fn parse_double_value(s: &str) -> Result<f64> {
    parse_double(s).ok_or_else(|| {
        if s.is_empty() { anyhow!("java.lang.NumberFormatException: empty String") } else { anyhow!("java.lang.NumberFormatException: For input string: \"{s}\"") }
    })
}

/// Files.readAllLines with the header's index map; Java fails on a file without a header line
pub fn read_table(path: &str) -> Result<(Vec<String>, HashMap<String, usize>)> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else {
        bail!("{path}: java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0")
    };
    let fields = fields_index_map(header, '\t');
    Ok((lines, fields))
}

// ---------------------------------------------------------------------------------------------
// LinxSvAnnotation

#[derive(Clone, Debug, PartialEq)]
pub struct LinxSvAnnotation {
    pub vcf_id_start: String,
    pub vcf_id_end: String,
    pub sv_id: i32,
    pub coords_start: String,
    pub coords_end: String,
    /// StructuralVariantType name
    pub sv_type: String,
    pub cluster_id: i32,
    pub cluster_reason: String,
    pub fragile_site_start: bool,
    pub fragile_site_end: bool,
    pub is_foldback: bool,
    pub line_type_start: String,
    pub line_type_end: String,
    pub junction_copy_number_min: f64,
    pub junction_copy_number_max: f64,
    pub gene_start: String,
    pub gene_end: String,
    pub local_topology_id_start: i32,
    pub local_topology_id_end: i32,
    pub local_topology_start: String,
    pub local_topology_end: String,
    pub local_ti_count_start: i32,
    pub local_ti_count_end: i32,
}

/// LinxSvAnnotation.generateFilename
pub fn sv_annotation_filename(base_path: &str, sample: &str, germline: bool) -> String {
    format!("{base_path}/{sample}{}", if germline { ".linx.germline.svs.tsv" } else { ".linx.svs.tsv" })
}

/// LinxSvAnnotation.read
pub fn read_sv_annotations(path: &str) -> Result<Vec<LinxSvAnnotation>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split(line, '\t') };
        let vcf_id_start = match f.get("vcfIdStart") {
            Some(&i) => r.at(i)?.to_string(),
            None => r.string("vcfId")?,
        };
        let vcf_id_end = match f.get("vcfIdEnd") {
            Some(&i) => r.at(i)?.to_string(),
            None => vcf_id_start.clone(),
        };
        let sv_type = match f.get("type") {
            Some(&i) => enum_value(r.at(i)?, &SV_TYPES, "StructuralVariantType")?,
            None => "BND".to_string(),
        };
        out.push(LinxSvAnnotation {
            vcf_id_start,
            vcf_id_end,
            sv_id: r.int("svId")?,
            coords_start: r.str_or("coordsStart", "")?.to_string(),
            coords_end: r.str_or("coordsEnd", "")?.to_string(),
            sv_type,
            cluster_id: r.int("clusterId")?,
            cluster_reason: r.string("clusterReason")?,
            fragile_site_start: r.bool("fragileSiteStart")?,
            fragile_site_end: r.bool("fragileSiteEnd")?,
            is_foldback: r.bool("isFoldback")?,
            line_type_start: r.string("lineTypeStart")?,
            line_type_end: r.string("lineTypeEnd")?,
            junction_copy_number_min: r.double("junctionCopyNumberMin")?,
            junction_copy_number_max: r.double("junctionCopyNumberMax")?,
            gene_start: r.string("geneStart")?,
            gene_end: r.string("geneEnd")?,
            local_topology_id_start: r.int("localTopologyIdStart")?,
            local_topology_id_end: r.int("localTopologyIdEnd")?,
            local_topology_start: r.string("localTopologyStart")?,
            local_topology_end: r.string("localTopologyEnd")?,
            local_ti_count_start: r.int("localTICountStart")?,
            local_ti_count_end: r.int("localTICountEnd")?,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// LinxBreakend

#[derive(Clone, Debug, PartialEq)]
pub struct LinxBreakend {
    pub id: i32,
    pub sv_id: i32,
    pub vcf_id: String,
    pub coords: String,
    pub is_start: bool,
    pub gene: String,
    pub transcript_id: String,
    pub canonical: bool,
    pub gene_orientation: String,
    pub disruptive: bool,
    /// ReportedStatus name
    pub reported_status: String,
    pub undisrupted_copy_number: f64,
    pub region_type: String,
    pub coding_type: String,
    pub biotype: String,
    pub exonic_base_phase: i32,
    pub next_splice_exon_rank: i32,
    pub next_splice_exon_phase: i32,
    pub next_splice_distance: i32,
    pub total_exon_count: i32,
    pub exon_up: i32,
    pub exon_down: i32,
}

/// LinxBreakend.generateFilename
pub fn breakend_filename(base_path: &str, sample: &str, germline: bool) -> String {
    format!("{base_path}/{sample}{}", if germline { ".linx.germline.breakend.tsv" } else { ".linx.breakend.tsv" })
}

/// LinxBreakend.read
pub fn read_breakends(path: &str) -> Result<Vec<LinxBreakend>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split(line, '\t') };
        let reported_status = match f.get("reportedStatus") {
            Some(&i) => enum_value(r.at(i)?, &REPORTED_STATUSES, "ReportedStatus")?,
            None => {
                if r.bool("reportedDisruption")? { "REPORTED".to_string() } else { "NONE".to_string() }
            }
        };
        out.push(LinxBreakend {
            id: r.int("id")?,
            sv_id: r.int("svId")?,
            vcf_id: r.str_or("vcfId", "")?.to_string(),
            coords: r.str_or("coords", "")?.to_string(),
            is_start: r.bool("isStart")?,
            gene: r.string("gene")?,
            transcript_id: r.string("transcriptId")?,
            canonical: r.bool("canonical")?,
            gene_orientation: r.string("geneOrientation")?,
            disruptive: r.bool("disruptive")?,
            reported_status,
            undisrupted_copy_number: r.double("undisruptedCopyNumber")?,
            region_type: r.enum_value("regionType", &TRANSCRIPT_REGION_TYPES, "TranscriptRegionType")?,
            coding_type: r.enum_value("codingType", &TRANSCRIPT_CODING_TYPES, "TranscriptCodingType")?,
            biotype: r.string("biotype")?,
            exonic_base_phase: r.int("exonicBasePhase")?,
            next_splice_exon_rank: r.int("nextSpliceExonRank")?,
            next_splice_exon_phase: r.int("nextSpliceExonPhase")?,
            next_splice_distance: r.int("nextSpliceDistance")?,
            total_exon_count: r.int("totalExonCount")?,
            exon_up: r.int("exonUp")?,
            exon_down: r.int("exonDown")?,
        });
    }
    Ok(out)
}

/// LinxBreakend.chromosomeFromCoords: `coords.split(":", 3)[0]`
pub fn chromosome_from_coords(coords: &str) -> String {
    coords.splitn(3, ':').next().unwrap_or("").to_string()
}

/// LinxBreakend.positionFromCoords
pub fn position_from_coords(coords: &str) -> Result<i32> {
    let parts: Vec<&str> = coords.splitn(3, ':').collect();
    let p = parts.get(1).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: coords({coords})"))?;
    parse_int_value(p)
}

/// LinxBreakend.orientationFromCoords
pub fn orientation_from_coords(coords: &str) -> Result<i8> {
    let parts: Vec<&str> = coords.splitn(3, ':').collect();
    let o = parts.get(2).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: coords({coords})"))?;
    parse_byte(o).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{o}\""))
}

// ---------------------------------------------------------------------------------------------
// LinxCluster

#[derive(Clone, Debug, PartialEq)]
pub struct LinxCluster {
    pub cluster_id: i32,
    pub category: String,
    pub synthetic: bool,
    pub resolved_type: String,
    pub cluster_count: i32,
    pub cluster_desc: String,
}

/// LinxCluster.generateFilename
pub fn cluster_filename(base_path: &str, sample: &str, germline: bool) -> String {
    format!("{base_path}/{sample}{}", if germline { ".linx.germline.clusters.tsv" } else { ".linx.clusters.tsv" })
}

/// LinxCluster.read
pub fn read_clusters(path: &str) -> Result<Vec<LinxCluster>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split(line, '\t') };
        out.push(LinxCluster {
            cluster_id: r.int("clusterId")?,
            category: r.string("category")?,
            synthetic: r.bool("synthetic")?,
            resolved_type: r.string("resolvedType")?,
            cluster_count: r.int("clusterCount")?,
            cluster_desc: r.string("clusterDesc")?,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// LinxFusion

#[derive(Clone, Debug, PartialEq)]
pub struct LinxFusion {
    pub five_prime_breakend_id: i32,
    pub three_prime_breakend_id: i32,
    pub name: String,
    pub reported: bool,
    /// KnownFusionType name as written (ORANGE converts it with LinxFusionType.valueOf)
    pub reported_type: String,
    pub reportable_reasons: Vec<String>,
    pub phased: String,
    pub likelihood: String,
    pub five_prime_vcf_id: String,
    pub three_prime_vcf_id: String,
    pub five_prime_coords: String,
    pub three_prime_coords: String,
    pub chain_length: i32,
    pub chain_links: i32,
    pub chain_terminated: bool,
    pub domains_kept: String,
    pub domains_lost: String,
    pub skipped_exons_up: i32,
    pub skipped_exons_down: i32,
    pub fused_exon_up: i32,
    pub fused_exon_down: i32,
}

/// LinxFusion.generateFilename
pub fn fusion_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.linx.fusion.tsv")
}

/// FileDelimiters.splitEnumsFromStr(str, FusionReportableReason): `str.split(";")`, each through valueOf
/// (so an empty string fails, as `"".split(";")` is `[""]`)
pub fn reportable_reasons_from_str(s: &str) -> Result<Vec<String>> {
    split(s, ';').into_iter().map(|n| enum_value(n, &FUSION_REPORTABLE_REASONS, "FusionReportableReason")).collect()
}

/// LinxFusion.read
pub fn read_fusions(path: &str) -> Result<Vec<LinxFusion>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split(line, '\t') };
        let five_prime_breakend_id = r.int("fivePrimeBreakendId")?;
        let three_prime_breakend_id = r.int("threePrimeBreakendId")?;
        let name = r.string("name")?;
        let reported = r.bool("reported")?;
        let reported_type = r.string("reportedType")?;
        let reportable_reasons = match f.get("reportableReasons") {
            Some(&i) => reportable_reasons_from_str(r.at(i)?)?,
            None => Vec::new(),
        };
        out.push(LinxFusion {
            five_prime_breakend_id,
            three_prime_breakend_id,
            name,
            reported,
            reported_type,
            reportable_reasons,
            phased: r.enum_value("phased", &FUSION_PHASED_TYPES, "FusionPhasedType")?,
            likelihood: r.enum_value("likelihood", &FUSION_LIKELIHOOD_TYPES, "FusionLikelihoodType")?,
            five_prime_vcf_id: r.str_or("fivePrimeVcfId", "")?.to_string(),
            three_prime_vcf_id: r.str_or("threePrimeVcfId", "")?.to_string(),
            five_prime_coords: r.str_or("fivePrimeCoords", "")?.to_string(),
            three_prime_coords: r.str_or("threePrimeCoords", "")?.to_string(),
            chain_length: r.int("chainLength")?,
            chain_links: r.int("chainLinks")?,
            chain_terminated: r.bool("chainTerminated")?,
            domains_kept: r.string("domainsKept")?,
            domains_lost: r.string("domainsLost")?,
            skipped_exons_up: r.int("skippedExonsUp")?,
            skipped_exons_down: r.int("skippedExonsDown")?,
            fused_exon_up: r.int("fusedExonUp")?,
            fused_exon_down: r.int("fusedExonDown")?,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// LinxDriver

#[derive(Clone, Debug, PartialEq)]
pub struct LinxDriver {
    pub cluster_id: i32,
    pub gene: String,
    /// DriverEventType name
    pub event_type: String,
}

/// LinxDriver.generateFilename
pub fn driver_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.linx.drivers.tsv")
}

/// LinxDriver.generateCatalogFilename
pub fn driver_catalog_filename(base_path: &str, sample: &str, somatic: bool) -> String {
    format!("{base_path}/{sample}{}", if somatic { ".linx.driver.catalog.tsv" } else { ".linx.germline.driver.catalog.tsv" })
}

/// LinxDriver.read
pub fn read_drivers(path: &str) -> Result<Vec<LinxDriver>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        out.push(LinxDriver {
            cluster_id: r.int("clusterId")?,
            gene: r.string("gene")?,
            event_type: r.enum_value("eventType", &DRIVER_EVENT_TYPES, "DriverEventType")?,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// LinxGermlineDisruption (the fields Java parses; ORANGE keeps only svId)

#[derive(Clone, Debug, PartialEq)]
pub struct LinxGermlineDisruption {
    pub sv_id: i32,
    pub chromosome_start: String,
    pub chromosome_end: String,
    pub position_start: i32,
    pub position_end: i32,
    pub sv_type: String,
    pub gene: String,
    pub cluster_id: i32,
}

/// LinxGermlineDisruption.generateFilename
pub fn germline_disruption_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.linx.germline.disruption.tsv")
}

/// LinxGermlineDisruption.read: every field parsed in Java's order, so that a malformed file fails as in Java
pub fn read_germline_disruptions(path: &str) -> Result<Vec<LinxGermlineDisruption>> {
    let (lines, f) = read_table(path)?;
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let r = Row { file: path, fields: &f, values: split_all(line, '\t') };
        let sv_id = r.int_or("svId", 0)?;
        r.str_or("vcfId", "")?;
        let chromosome_start = r.string("chromosomeStart")?;
        let chromosome_end = r.string("chromosomeEnd")?;
        let position_start = r.int_or("positionStart", 0)?;
        let position_end = r.int_or("positionEnd", 0)?;
        r.byte("orientStart")?;
        r.byte("orientEnd")?;
        let sv_type = r.enum_value("type", &SV_TYPES, "StructuralVariantType")?;
        r.str("filter")?;
        r.str("event")?;
        r.double_or("qualScore", 0.0)?;
        r.str_or("homologySequenceStart", "")?;
        r.str_or("homologySequenceEnd", "")?;
        r.int_or("inexactHomologyStart", 0)?;
        r.int_or("inexactHomologyEnd", 0)?;
        for name in [
            "junctionCopyNumber",
            "adjustedAFStart",
            "adjustedAFEnd",
            "adjustedCopyNumberStart",
            "adjustedCopyNumberEnd",
            "adjustedCopyNumberChangeStart",
            "adjustedCopyNumberChangeEnd",
        ] {
            r.double_or(name, 0.0)?;
        }
        for name in [
            "germlineFragments",
            "germlineReferenceFragmentsStart",
            "germlineReferenceFragmentsEnd",
            "tumorFragments",
            "tumorReferenceFragmentsStart",
            "tumorReferenceFragmentsEnd",
        ] {
            r.int_or(name, 0)?;
        }
        for name in ["insertSequence", "insertSequenceAlignments", "insertSequenceRepeatClass", "insertSequenceRepeatType"] {
            r.str(name)?;
        }
        let gene = r.str_or("gene", "")?.to_string();
        let cluster_id = r.int_or("clusterId", 0)?;
        r.int_or("clusterCount", 0)?;
        r.str("resolvedType")?;
        r.str("linkedByStart")?;
        r.str("linkedByEnd")?;
        r.int_or("cohortFrequency", 0)?;
        out.push(LinxGermlineDisruption {
            sv_id,
            chromosome_start,
            chromosome_end,
            position_start,
            position_end,
            sv_type,
            gene,
            cluster_id,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &std::path::Path, name: &str, text: &str) -> String {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn coords() {
        assert_eq!(chromosome_from_coords("chr9:21972347:-1"), "chr9");
        assert_eq!(position_from_coords("chr9:21972347:-1").unwrap(), 21972347);
        assert_eq!(orientation_from_coords("chr9:21972347:-1").unwrap(), -1);
        assert!(position_from_coords("chr9").is_err());
        assert_eq!(vis_plot_cluster_prefix("12"), ".cluster-12");
    }

    #[test]
    fn readers() {
        let dir = std::env::temp_dir().join(format!("linx_files_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = write(
            &dir,
            "S.linx.breakend.tsv",
            "id\tsvId\tvcfId\tcoords\tisStart\tgene\ttranscriptId\tcanonical\tgeneOrientation\tdisruptive\treportedStatus\tundisruptedCopyNumber\tregionType\tcodingType\tbiotype\texonicBasePhase\tnextSpliceExonRank\tnextSpliceExonPhase\tnextSpliceDistance\ttotalExonCount\texonUp\texonDown\n\
             1\t2\tv\tchr1:100:1\ttrue\tGENE\tENST1\ttrue\tUpstream\ttrue\tREPORTED\t0.5\tINTRONIC\tCODING\tprotein_coding\t-1\t2\t0\t100\t5\t1\t2\n",
        );
        let b = read_breakends(&p).unwrap();
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].reported_status, "REPORTED");
        assert_eq!(b[0].exon_down, 2);
        // a trailing empty field is dropped by String.split, so the last column is out of range
        let p = write(&dir, "S.linx.clusters.tsv", "clusterId\tcategory\tsynthetic\tresolvedType\tclusterCount\tclusterDesc\n3\tSIMPLE\tfalse\tDEL\t1\t\n");
        assert!(read_clusters(&p).is_err());
        let p = write(&dir, "S.linx.drivers.tsv", "clusterId\tgene\teventType\n3\tG\tDEL\n4\tH\tNOPE\n");
        assert!(read_drivers(&p).is_err());
        assert!(reportable_reasons_from_str("").is_err());
        assert_eq!(reportable_reasons_from_str("OK;NMD").unwrap(), vec!["OK", "NMD"]);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
