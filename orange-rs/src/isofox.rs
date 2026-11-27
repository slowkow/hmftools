//! ISOFOX's part of the report: IsofoxDataLoader, IsofoxInterpreter and IsofoxConversion (orange/algo/isofox).

use anyhow::{bail, Result};
use hmf_common::rna;

use crate::model::{self, enum_set, LinxFusion, RNA_QC_STATUSES};

pub struct IsofoxData {
    pub summary: rna::RnaStatistics,
    pub gene_expressions: Vec<rna::GeneExpression>,
    pub fusions: Vec<rna::RnaFusion>,
    pub novel_splice_junctions: Vec<rna::NovelSpliceJunction>,
}

fn mandatory(path: String) -> Result<String> {
    if !std::path::Path::new(&path).exists() {
        bail!("java.lang.IllegalArgumentException: Path [{path}] did not exist and is mandatory.");
    }
    Ok(path)
}

/// IsofoxDataLoader.checkRnaVsTumorId: the RNA id when its summary exists, else the tumor id when that one's
/// does (backwards compatibility), else None
pub fn check_rna_vs_tumor_id(rna_id: Option<&str>, tumor_id: &str, dir: &str) -> Option<String> {
    let exists = |id: &str| std::path::Path::new(&rna::rna_statistics_filename(dir, id)).exists();
    let rna_name = rna_id.unwrap_or("null");
    if exists(rna_name) {
        return Some(rna_name.to_string());
    }
    if exists(tumor_id) {
        return Some(tumor_id.to_string());
    }
    None
}

/// IsofoxDataLoader.load
pub fn load(rna_id: &str, dir: &str) -> Result<IsofoxData> {
    let summary = rna::read_rna_statistics(&mandatory(rna::rna_statistics_filename(dir, rna_id))?)?;
    let gene_expressions = rna::read_gene_expressions(&mandatory(rna::gene_expression_filename(dir, rna_id))?)?;
    let fusions = rna::read_rna_fusions(&mandatory(rna::rna_fusion_filename(dir, rna_id))?)?;
    let novel_splice_junctions = rna::read_novel_splice_junctions(&mandatory(rna::novel_splice_junction_filename(dir, rna_id))?)?;
    Ok(IsofoxData { summary, gene_expressions, fusions, novel_splice_junctions })
}

/// IsofoxConversion.convert(GeneExpression)
pub fn convert_gene_expression(e: &rna::GeneExpression) -> model::GeneExpression {
    let cancer = e.median_tpm_cancer > 0.0;
    model::GeneExpression {
        gene: e.gene_name.clone(),
        tpm: e.tpm,
        median_tpm_cohort: e.median_tpm_cohort,
        percentile_cohort: e.percentile_cohort,
        median_tpm_cancer: if cancer { Some(e.median_tpm_cancer) } else { None },
        percentile_cancer: if cancer { Some(e.percentile_cancer) } else { None },
    }
}

/// LinxFusion.geneNames: `name.split("_", 2)`
fn gene_names(name: &str) -> Vec<&str> {
    name.splitn(2, '_').collect()
}

/// IsofoxConversion.convert(RnaFusion)
pub(crate) fn convert_fusion(f: &rna::RnaFusion) -> model::RnaFusion {
    let names = gene_names(&f.name);
    let known_type = if names.len() == 2 && names[0] == names[1] {
        "EXON_DEL_DUP"
    } else {
        match f.known_type.as_str() {
            "KNOWN_PAIR" => "KNOWN_PAIR",
            "PROMISCUOUS_5" => "PROMISCUOUS_5",
            "PROMISCUOUS_3" => "PROMISCUOUS_3",
            _ => "NONE",
        }
    };
    model::RnaFusion {
        gene_start: names.first().map(|s| s.to_string()),
        gene_end: names.get(1).map(|s| s.to_string()),
        chromosome_start: f.chromosome_up.clone(),
        chromosome_end: f.chromosome_down.clone(),
        position_start: f.position_up,
        position_end: f.position_down,
        junction_type_start: f.junction_type_up.clone(),
        junction_type_end: f.junction_type_down.clone(),
        known_type: known_type.to_string(),
        sv_type: f.sv_type.clone(),
        split_fragments: f.split_fragments,
        realigned_frags: f.realigned_frags,
        discordant_frags: f.discordant_frags,
        depth_start: f.depth_up,
        depth_end: f.depth_down,
        cohort_frequency: f.cohort_frequency,
    }
}

/// IsofoxInterpreter.findFusions: known pairs that LINX did not find in the DNA
pub fn find_fusions(rna_fusions: &[rna::RnaFusion], linx_fusions: &[LinxFusion]) -> Vec<model::RnaFusion> {
    let mut out = Vec::new();
    for f in rna_fusions {
        let names = gene_names(&f.name);
        if names.len() != 2 || names[0].is_empty() || names[1].is_empty() {
            continue;
        }
        if linx_fusions.iter().any(|l| l.gene_up == names[0] && l.gene_down == names[1]) {
            continue;
        }
        if f.known_type == "KNOWN_PAIR" {
            out.push(convert_fusion(f));
        }
    }
    out
}

/// IsofoxConversion.convert(NovelSpliceJunction): AltSpliceJunctionType.valueOf (the data model has no CANONICAL)
fn convert_junction(j: &rna::NovelSpliceJunction) -> Result<model::NovelSpliceJunction> {
    if j.junction_type == "CANONICAL" {
        bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType.CANONICAL");
    }
    Ok(model::NovelSpliceJunction {
        gene: j.gene_name.clone(),
        chromosome: j.chromosome.clone(),
        junction_start: j.junction_start,
        junction_end: j.junction_end,
        junction_type: j.junction_type.clone(),
        exon_start: j.exon_start,
        exon_end: j.exon_end,
        fragment_count: j.fragment_count,
        depth_start: j.depth_start,
        depth_end: j.depth_end,
        region_start: j.region_start.clone(),
        region_end: j.region_end.clone(),
        cohort_frequency: j.cohort_frequency,
    })
}

/// IsofoxInterpreter.interpret
pub fn interpret(data: &IsofoxData, linx_fusions: &[LinxFusion]) -> Result<model::IsofoxRecord> {
    let mut high = Vec::new();
    // findExpressionOutliers: the reported genes with cohort data; the low list is never filled
    for e in &data.gene_expressions {
        if e.median_tpm_cohort <= 0.0 {
            continue;
        }
        if e.reported_status == "REPORTED" {
            high.push(convert_gene_expression(e));
        }
    }
    let fusions = find_fusions(&data.fusions, linx_fusions);
    let junctions = data.novel_splice_junctions.iter().map(convert_junction).collect::<Result<Vec<_>>>()?;
    let s = &data.summary;
    Ok(model::IsofoxRecord {
        summary: model::RnaStatistics {
            qc_status: enum_set(&s.qc_status, &RNA_QC_STATUSES),
            total_fragments: s.total_fragments,
            duplicate_fragments: s.duplicate_fragments,
            spliced_fragment_perc: s.spliced_fragment_perc,
            unspliced_fragment_perc: s.unspliced_fragment_perc,
            alt_fragment_perc: s.alt_fragment_perc,
            chimeric_fragment_perc: s.chimeric_fragment_perc,
        },
        high_expression_genes: high,
        low_expression_genes: Vec::new(),
        fusions,
        novel_splice_junctions: junctions,
    })
}
