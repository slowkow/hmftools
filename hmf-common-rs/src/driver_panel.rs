//! hmf-common `driver.panel`: the driver gene panel file (DriverGeneFile, DriverGene, DriverGenePanelConfig)
//! and the regions of the panel's genes (DriverGeneRegions.buildDriverGeneRegions and
//! getTranscriptRegions).

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashSet;

use crate::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use crate::java_hash::JavaHashMap;
use crate::java_text::{self, parse_bool, parse_double};
use crate::regions::positions_overlap;
use crate::segmentation::HumanChromosome;

pub const DRIVER_GENE_PANEL: &str = "driver_gene_panel";
pub const DRIVER_GENE_PANEL_DESC: &str = "Path to driver gene panel";

/// PurpleCommon defaults for the optional columns
pub const DEFAULT_DRIVER_HET_DELETION_THRESHOLD: f64 = 0.6;
pub const DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO: f64 = 3.0;

/// DriverCategory
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DriverCategory {
    Onco,
    Tsg,
}

impl DriverCategory {
    /// valueOf: the exact name
    pub fn value_of(s: &str) -> Result<Self> {
        match s {
            "ONCO" => Ok(DriverCategory::Onco),
            "TSG" => Ok(DriverCategory::Tsg),
            _ => bail!("No enum constant DriverCategory.{s}"),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            DriverCategory::Onco => "ONCO",
            DriverCategory::Tsg => "TSG",
        }
    }
}

/// DriverGeneGermlineReporting
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GermlineReporting {
    Any,
    VariantNotLost,
    WildtypeLost,
    None,
}

impl GermlineReporting {
    /// valueOf: the exact name
    pub fn value_of(s: &str) -> Result<Self> {
        match s {
            "ANY" => Ok(GermlineReporting::Any),
            "VARIANT_NOT_LOST" => Ok(GermlineReporting::VariantNotLost),
            "WILDTYPE_LOST" => Ok(GermlineReporting::WildtypeLost),
            "NONE" => Ok(GermlineReporting::None),
            _ => bail!("No enum constant DriverGeneGermlineReporting.{s}"),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            GermlineReporting::Any => "ANY",
            GermlineReporting::VariantNotLost => "VARIANT_NOT_LOST",
            GermlineReporting::WildtypeLost => "WILDTYPE_LOST",
            GermlineReporting::None => "NONE",
        }
    }
}

/// DriverGene
#[derive(Clone, Debug, PartialEq)]
pub struct DriverGene {
    pub gene: String,
    pub report_missense_and_inframe: bool,
    pub report_nonsense_and_frameshift: bool,
    pub report_splice: bool,
    pub report_deletion: bool,
    pub report_het_deletion: bool,
    pub report_loh: bool,
    pub het_deletion_threshold: f64,
    pub report_disruption: bool,
    pub report_amplification: bool,
    pub amplification_ratio: f64,
    pub report_somatic_hotspot: bool,
    pub likelihood_type: DriverCategory,
    pub report_germline_variant: GermlineReporting,
    pub report_germline_hotspot: GermlineReporting,
    pub report_germline_disruption: GermlineReporting,
    pub report_germline_deletion: GermlineReporting,
    pub report_germline_amplification: bool,
    pub additional_reported_transcripts: Vec<String>,
    pub report_pgx: bool,
    pub report_high_expression: bool,
    pub report_low_expression: bool,
    pub report_novel_splice_junction: bool,
}

impl DriverGene {
    pub fn report_somatic(&self) -> bool {
        self.report_missense_and_inframe || self.report_nonsense_and_frameshift || self.report_splice || self.report_somatic_hotspot
    }
    pub fn report_germline(&self) -> bool {
        self.report_germline_variant != GermlineReporting::None
            || self.report_germline_hotspot != GermlineReporting::None
            || self.report_germline_deletion != GermlineReporting::None
            || self.report_germline_disruption != GermlineReporting::None
    }
}

/// The columns DriverGeneFile writes, in order
const COLUMNS: [&str; 23] = [
    "gene",
    "reportMissense",
    "reportNonsense",
    "reportSplice",
    "reportDeletion",
    "reportHetDeletion",
    "reportLoh",
    "hetDeletionThreshold",
    "reportDisruption",
    "reportAmplification",
    "amplificationRatio",
    "reportSomaticHotspot",
    "likelihoodType",
    "reportGermlineVariant",
    "reportGermlineHotspot",
    "reportGermlineDisruption",
    "reportGermlineDeletion",
    "reportGermlineAmplification",
    "additionalReportedTranscripts",
    "reportPGX",
    "reportHighExpression",
    "reportLowExpression",
    "reportNovelSpliceJunction",
];

/// DriverGeneFile.read(filename)
pub fn read_driver_genes(filename: &str) -> Result<Vec<DriverGene>> {
    let lines = java_text::read_lines(filename).with_context(|| format!("reading driver gene panel {filename}"))?;
    from_lines(&lines).with_context(|| format!("parsing driver gene panel {filename}"))
}

/// DriverGenePanelConfig.loadDriverGenes: the panel's genes, an empty list when the file cannot be read
/// (Java catches only the IOException; a malformed line is an error)
pub fn load_driver_genes(filename: Option<&str>) -> Result<Vec<DriverGene>> {
    let Some(f) = filename else { return Ok(Vec::new()) };
    let lines = match java_text::read_lines(f) {
        Ok(l) => l,
        Err(_) => return Ok(Vec::new()),
    };
    from_lines(&lines).with_context(|| format!("parsing driver gene panel {f}"))
}

/// DriverGeneFile.fromLines: the header names the columns; reportHetDeletion, reportLoh,
/// hetDeletionThreshold, amplificationRatio, reportGermlineAmplification and the RNA columns may be
/// absent
pub fn from_lines(lines: &[String]) -> Result<Vec<DriverGene>> {
    let header = lines.first().ok_or_else(|| anyhow!("empty driver gene panel"))?;
    let f = java_text::fields_index_map(header, '\t');
    let req = |name: &str| f.get(name).copied().ok_or_else(|| anyhow!("missing column {name}"));
    let opt = |name: &str| f.get(name).copied();
    let gene_i = req("gene")?;
    let missense_i = req("reportMissense")?;
    let nonsense_i = req("reportNonsense")?;
    let splice_i = req("reportSplice")?;
    let deletion_i = req("reportDeletion")?;
    let het_del_i = opt("reportHetDeletion");
    let loh_i = opt("reportLoh");
    let het_del_threshold_i = opt("hetDeletionThreshold");
    let disruption_i = req("reportDisruption")?;
    let amp_i = req("reportAmplification")?;
    let amp_ratio_i = opt("amplificationRatio");
    let hotspot_i = req("reportSomaticHotspot")?;
    let likelihood_i = req("likelihoodType")?;
    let germline_variant_i = req("reportGermlineVariant")?;
    let germline_hotspot_i = req("reportGermlineHotspot")?;
    let germline_disruption_i = req("reportGermlineDisruption")?;
    let germline_deletion_i = req("reportGermlineDeletion")?;
    let germline_amp_i = opt("reportGermlineAmplification");
    let alt_trans_i = req("additionalReportedTranscripts")?;
    let pgx_i = req("reportPGX")?;
    let low_exp_i = opt("reportLowExpression");
    let high_exp_i = opt("reportHighExpression");
    let novel_sj_i = opt("reportNovelSpliceJunction");

    let mut genes = Vec::with_capacity(lines.len().saturating_sub(1));
    for line in &lines[1..] {
        let v = java_text::split_all(line, '\t');
        let val = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("ArrayIndexOutOfBounds in line: {line}"));
        let opt_bool = |i: Option<usize>| -> Result<bool> { Ok(match i { Some(i) => parse_bool(val(i)?), None => false }) };
        let opt_double = |i: Option<usize>, default: f64| -> Result<f64> {
            Ok(match i {
                Some(i) => parse_double(val(i)?).ok_or_else(|| anyhow!("invalid number in line: {line}"))?,
                None => default,
            })
        };
        let other_trans: Vec<String> = java_text::split(val(alt_trans_i)?, ';').into_iter().filter(|x| !x.is_empty()).map(|x| x.to_string()).collect();
        // backwards compatibility prior to pipeline v5.32: true / false for the germline disruption
        let disruption_str = val(germline_disruption_i)?;
        let germline_disruption = match disruption_str.to_lowercase().as_str() {
            "true" => GermlineReporting::Any,
            "false" => GermlineReporting::None,
            _ => GermlineReporting::value_of(disruption_str)?,
        };
        genes.push(DriverGene {
            gene: val(gene_i)?.to_string(),
            report_missense_and_inframe: parse_bool(val(missense_i)?),
            report_nonsense_and_frameshift: parse_bool(val(nonsense_i)?),
            report_splice: parse_bool(val(splice_i)?),
            report_deletion: parse_bool(val(deletion_i)?),
            report_het_deletion: opt_bool(het_del_i)?,
            report_loh: opt_bool(loh_i)?,
            het_deletion_threshold: opt_double(het_del_threshold_i, DEFAULT_DRIVER_HET_DELETION_THRESHOLD)?,
            report_disruption: parse_bool(val(disruption_i)?),
            report_amplification: parse_bool(val(amp_i)?),
            amplification_ratio: opt_double(amp_ratio_i, DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)?,
            report_somatic_hotspot: parse_bool(val(hotspot_i)?),
            likelihood_type: DriverCategory::value_of(val(likelihood_i)?)?,
            report_germline_variant: GermlineReporting::value_of(&val(germline_variant_i)?.to_uppercase())?,
            report_germline_hotspot: GermlineReporting::value_of(&val(germline_hotspot_i)?.to_uppercase())?,
            report_germline_disruption: germline_disruption,
            report_germline_deletion: GermlineReporting::value_of(val(germline_deletion_i)?)?,
            report_germline_amplification: opt_bool(germline_amp_i)?,
            additional_reported_transcripts: other_trans,
            report_pgx: parse_bool(val(pgx_i)?),
            report_low_expression: opt_bool(low_exp_i)?,
            report_high_expression: opt_bool(high_exp_i)?,
            report_novel_splice_junction: opt_bool(novel_sj_i)?,
        });
    }
    Ok(genes)
}

/// DriverGeneFile.write with each value under its own column: the 23 columns, booleans as
/// "true"/"false", decimals as Double.toString.
///
/// Java's writer (`write_driver_genes_java`) puts reportLowExpression's value under the
/// reportHighExpression header and the other way round (its header follows the Columns enum, its values
/// the order of DriverGeneFile.toString), so a panel it writes reads back with the two swapped.
pub fn write_driver_genes(filename: &str, genes: &[DriverGene]) -> Result<()> {
    write_genes(filename, genes, false)
}

/// DriverGeneFile.write as Java has it: the low and high expression values swapped (see `write_driver_genes`)
pub fn write_driver_genes_java(filename: &str, genes: &[DriverGene]) -> Result<()> {
    write_genes(filename, genes, true)
}

fn write_genes(filename: &str, genes: &[DriverGene], java_swap: bool) -> Result<()> {
    let mut out = COLUMNS.join("\t");
    out.push('\n');
    let b = |x: bool| if x { "true" } else { "false" };
    for g in genes {
        let fields: Vec<String> = vec![
            g.gene.clone(),
            b(g.report_missense_and_inframe).into(),
            b(g.report_nonsense_and_frameshift).into(),
            b(g.report_splice).into(),
            b(g.report_deletion).into(),
            b(g.report_het_deletion).into(),
            b(g.report_loh).into(),
            crate::fmt::java_double_to_string(g.het_deletion_threshold),
            b(g.report_disruption).into(),
            b(g.report_amplification).into(),
            crate::fmt::java_double_to_string(g.amplification_ratio),
            b(g.report_somatic_hotspot).into(),
            g.likelihood_type.name().into(),
            g.report_germline_variant.name().into(),
            g.report_germline_hotspot.name().into(),
            g.report_germline_disruption.name().into(),
            g.report_germline_deletion.name().into(),
            b(g.report_germline_amplification).into(),
            g.additional_reported_transcripts.join(";"),
            b(g.report_pgx).into(),
            b(if java_swap { g.report_low_expression } else { g.report_high_expression }).into(),
            b(if java_swap { g.report_high_expression } else { g.report_low_expression }).into(),
            b(g.report_novel_splice_junction).into(),
        ];
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    std::fs::write(filename, out).with_context(|| format!("writing {filename}"))
}

/// GeneRegion: a gene's exon region with its name and exon rank
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub gene_name: String,
    pub exon_rank: i32,
}

impl GeneRegion {
    pub fn base_length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// DriverGeneRegions.findInvalidDriverGenes: the genes the cache has no gene of that name for
pub fn find_invalid_driver_genes<'a>(genes: &'a [DriverGene], cache: &EnsemblDataCache) -> Vec<&'a DriverGene> {
    genes.iter().filter(|g| cache.gene_data_by_name(&g.gene).is_none()).collect()
}

const SPLICE_SIZE: i32 = 10;

/// DriverGeneRegions.getTranscriptRegions: each exon (widened by 10 bases on the inner sides when
/// `include_splice`) within the coding region (or the whole transcript when `include_utr` or non-coding)
pub fn transcript_regions(gene: &GeneData, trans: &TranscriptData, include_utr: bool, include_splice: bool) -> Vec<GeneRegion> {
    let whole = include_utr || trans.non_coding();
    let start_position = if whole { trans.trans_start } else { trans.coding_start.unwrap() };
    let end_position = if whole { trans.trans_end } else { trans.coding_end.unwrap_or(trans.trans_end) };
    let n = trans.exons.len();
    let mut regions = Vec::new();
    for (i, exon) in trans.exons.iter().enumerate() {
        let (mut exon_start, mut exon_end) = (exon.start, exon.end);
        if include_splice {
            if i > 0 {
                exon_start -= SPLICE_SIZE;
            }
            if i + 1 < n {
                exon_end += SPLICE_SIZE;
            }
        }
        if positions_overlap(start_position, end_position, exon_start, exon_end) {
            regions.push(GeneRegion {
                chromosome: gene.chromosome.clone(),
                start: start_position.max(exon_start),
                end: end_position.min(exon_end),
                gene_name: gene.gene_name.clone(),
                exon_rank: exon.rank,
            });
        }
    }
    regions
}

/// DriverGeneRegions.buildDriverGeneRegions: for each human chromosome (1-22, X, Y, versioned), the coding
/// exon regions with splice margins of the canonical transcript of every gene of the chromosome named
/// in `driver_genes`, sorted by start (stable); the map is put in chromosome order and iterates as Java's
/// HashMap. Java fails with a NullPointerException where a chromosome has no genes in the cache or a gene
/// has no canonical transcript; those are errors here.
pub fn build_driver_gene_regions(cache: &EnsemblDataCache, driver_genes: &[String]) -> Result<JavaHashMap<String, Vec<GeneRegion>>> {
    let gene_set: HashSet<&str> = driver_genes.iter().map(|s| s.as_str()).collect();
    let mut map: JavaHashMap<String, Vec<GeneRegion>> = JavaHashMap::new();
    for c in HumanChromosome::ALL {
        let chr = cache.ref_genome_version().versioned_chromosome(c.short_name());
        let genes = cache.chr_gene_data_map().get(&chr).ok_or_else(|| anyhow!("NullPointerException: no Ensembl genes for chromosome {chr}"))?;
        let mut regions: Vec<GeneRegion> = Vec::new();
        for g in genes {
            if !gene_set.contains(g.gene_name.as_str()) {
                continue;
            }
            let t = cache
                .canonical_transcript_data(&g.gene_id)
                .ok_or_else(|| anyhow!("NullPointerException: no canonical transcript for gene {}", g.gene_name))?;
            if !t.non_coding() && t.coding_end.is_none() {
                bail!("NullPointerException: no coding end for transcript {}", t.trans_name);
            }
            regions.extend(transcript_regions(g, t, false, true));
        }
        // Collections.sort by ChrBaseRegion.compareTo: one chromosome, so by start, stable
        regions.sort_by_key(|r| r.start);
        map.insert(chr, regions);
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::regions::RefGenomeVersion;

    const PANEL: &str = "gene\treportMissense\treportNonsense\treportSplice\treportDeletion\treportHetDeletion\thetDeletionThreshold\treportLoh\treportDisruption\treportAmplification\tamplificationRatio\treportSomaticHotspot\tlikelihoodType\treportGermlineVariant\treportGermlineHotspot\treportGermlineDisruption\treportGermlineDeletion\treportGermlineAmplification\tadditionalReportedTranscripts\treportPGX\treportHighExpression\treportNovelSpliceJunction\n\
        GENEA\tTRUE\tTRUE\tFALSE\tTRUE\tFALSE\t0.6\tFALSE\tTRUE\tFALSE\t3\tTRUE\tTSG\tany\tnone\ttrue\tNONE\tFALSE\tENST1;ENST2\tFALSE\tTRUE\tFALSE\n\
        GENEB\tFALSE\tFALSE\tFALSE\tFALSE\tTRUE\t0.5\tTRUE\tFALSE\tTRUE\t2.5\tFALSE\tONCO\tVARIANT_NOT_LOST\tWILDTYPE_LOST\tWILDTYPE_LOST\tANY\tTRUE\t\tTRUE\tFALSE\tTRUE\n";

    #[test]
    fn reads_and_writes_the_panel() {
        let lines: Vec<String> = PANEL.lines().map(|s| s.to_string()).collect();
        let genes = from_lines(&lines).unwrap();
        assert_eq!(genes.len(), 2);
        let a = &genes[0];
        assert!(a.report_missense_and_inframe && !a.report_splice && a.report_disruption);
        assert_eq!(a.report_germline_variant, GermlineReporting::Any);
        assert_eq!(a.report_germline_disruption, GermlineReporting::Any);
        assert_eq!(a.additional_reported_transcripts, vec!["ENST1", "ENST2"]);
        assert!(a.report_high_expression && !a.report_low_expression);
        let b = &genes[1];
        assert_eq!(b.likelihood_type, DriverCategory::Onco);
        assert_eq!(b.amplification_ratio, 2.5);
        assert!(b.report_loh && b.report_germline_amplification && b.report_novel_splice_junction);
        assert!(b.additional_reported_transcripts.is_empty());
        let dir = std::env::temp_dir().join(format!("hmf_common_panel_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("panel.tsv");
        write_driver_genes(p.to_str().unwrap(), &genes).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        assert!(text.lines().nth(1).unwrap().starts_with("GENEA\ttrue\ttrue\tfalse\ttrue\tfalse\tfalse\t0.6\ttrue\tfalse\t3.0\ttrue\tTSG\tANY\tNONE\tANY\tNONE\tfalse\tENST1;ENST2"));
        assert_eq!(read_driver_genes(p.to_str().unwrap()).unwrap(), genes);
        // Java's writer swaps the expression columns
        write_driver_genes_java(p.to_str().unwrap(), &genes).unwrap();
        let back = read_driver_genes(p.to_str().unwrap()).unwrap();
        assert!(!back[0].report_high_expression && back[0].report_low_expression);
        assert!(load_driver_genes(Some("/nonexistent/panel.tsv")).unwrap().is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn gene_regions() {
        let dir = std::env::temp_dir().join(format!("hmf_common_panel_regions_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut genes = String::from("GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,KaryotypeBand,Synonyms\n");
        for c in HumanChromosome::ALL {
            genes.push_str(&format!("G{0},GENE{0},{0},1,100,1000,p1,\n", c.short_name()));
        }
        std::fs::write(dir.join("ensembl_gene_data.38.csv"), genes).unwrap();
        // GENE1: three exons, coding 150-850; GENE2 non-coding
        let trans = "GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd,RefSeqId\n\
            G1,1,1,1,ENST1,protein_coding,100,900,1,100,200,-1,0,150,850,NULL\n\
            G1,1,1,1,ENST1,protein_coding,100,900,2,400,500,0,0,150,850,NULL\n\
            G1,1,1,1,ENST1,protein_coding,100,900,3,845,900,0,-1,150,850,NULL\n\
            G2,2,1,2,ENST2,lincRNA,100,900,1,100,300,-1,-1,NULL,NULL,NULL\n";
        std::fs::write(dir.join("ensembl_trans_exon_data.38.csv"), trans).unwrap();
        let mut cache = EnsemblDataCache::new(dir.to_str().unwrap(), RefGenomeVersion::V38);
        cache.set_required_data(true, false, false, true);
        cache.load(false).unwrap();
        cache.create_gene_name_id_map();
        let panel = vec!["GENE1".to_string(), "GENE2".to_string()];
        let m = build_driver_gene_regions(&cache, &panel).unwrap();
        assert_eq!(m.len(), 24);
        let r1 = m.get(&"chr1".to_string()).unwrap();
        assert_eq!(r1.iter().map(|r| (r.start, r.end, r.exon_rank)).collect::<Vec<_>>(), vec![(150, 210), (390, 510), (835, 850)].into_iter().zip([1, 2, 3]).map(|((a, b), c)| (a, b, c)).collect::<Vec<_>>());
        let r2 = m.get(&"chr2".to_string()).unwrap();
        assert_eq!((r2[0].start, r2[0].end), (100, 300));
        assert!(m.get(&"chrX".to_string()).unwrap().is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }
}
