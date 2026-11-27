//! GeneCopyNumber, GeneCopyNumberBuilder (the copy numbers over each transcript's exons) and
//! GeneCopyNumberFile.

use hmf_common::ensembl::{EnsemblDataCache, ExonData, GeneData, TranscriptData};

use crate::copynumber::PurpleCopyNumber;
use crate::types::*;

#[derive(Clone, Debug, PartialEq)]
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
    pub relative_min_copy_number: f64,
    pub somatic_regions: i32,
    pub min_regions: i32,
    pub min_region_start: i32,
    pub min_region_end: i32,
    pub depth_window_count: i32,
    pub gc_content: f64,
    pub min_region_start_support: SegmentSupport,
    pub min_region_end_support: SegmentSupport,
    pub min_region_method: CopyNumberMethod,
    pub driver_type: DriverType,
    pub reported_status: ReportedStatus,
}

impl GeneCopyNumber {
    pub fn min_region_bases(&self) -> i32 {
        self.min_region_end - self.min_region_start + 1
    }
}

fn weighted_average(cw: i64, cv: f64, nw: i64, nv: f64) -> f64 {
    if is_zero(cv) {
        return nv;
    }
    (cw as f64 * cv + nw as f64 * nv) / (cw + nw) as f64
}

struct Builder<'a> {
    copy_numbers: &'a [&'a PurpleCopyNumber],
    min_cn: f64,
    min_minor: f64,
    max_cn: f64,
    somatic_count: i32,
    previous: Option<&'a PurpleCopyNumber>,
    previous_cn: f64,
    exon: Option<&'a ExonData>,
    copy_number: Option<&'a PurpleCopyNumber>,
    min_regions: i32,
    min_region_start: i32,
    min_region_end: i32,
    min_start_support: SegmentSupport,
    min_end_support: SegmentSupport,
    min_method: CopyNumberMethod,
    depth_window_count: i32,
    min_gc: f64,
}

impl<'a> Builder<'a> {
    fn add_overlap(&mut self, exon: &ExonData, cn: &'a PurpleCopyNumber) {
        let min_end = exon.end.min(cn.end);
        let max_start = exon.start.max(cn.start);
        let overlap = 0.max(1 + min_end - max_start);
        if overlap <= 0 {
            return;
        }
        let current = cn.average_tumor_copy_number;
        self.max_cn = jmax(self.max_cn, current);
        self.min_minor = jmin(self.min_minor, cn.minor_allele_copy_number());
        if !equal(current, self.previous_cn) && !matches!(cn.method, CopyNumberMethod::GermlineHomDeletion | CopyNumberMethod::GermlineHet2HomDeletion) {
            self.somatic_count += 1;
        }
        let unprocessed = self.previous.is_none_or(|p| p != cn);
        if unprocessed {
            if less_than(current, self.min_cn) {
                self.min_regions = 1;
                self.min_cn = current;
                self.min_region_start = cn.start;
                self.min_start_support = cn.segment_start_support;
                self.min_region_end = cn.end;
                self.min_end_support = cn.segment_end_support;
                self.min_method = cn.method;
                self.depth_window_count = cn.depth_window_count;
                self.min_gc = cn.gc_content;
            } else if equal(current, self.min_cn) {
                self.min_region_end = cn.end;
                self.min_end_support = cn.segment_end_support;
                self.min_method = cn.method;
                if !equal(current, self.previous_cn) {
                    self.min_regions += 1;
                }
                self.min_gc = weighted_average(self.depth_window_count as i64, self.min_gc, cn.depth_window_count as i64, cn.gc_content);
                self.depth_window_count += cn.depth_window_count;
            }
        }
        self.previous_cn = current;
        self.previous = Some(cn);
    }

    fn create(mut self, gene: &GeneData, trans: &'a TranscriptData, ploidy: f64) -> GeneCopyNumber {
        let (mut ci, mut ei) = (0usize, 0usize);
        while ci < self.copy_numbers.len() || ei < trans.exons.len() {
            let cn = self.copy_numbers.get(ci).copied();
            let exon = trans.exons.get(ei);
            if cn.is_none() || exon.is_some_and(|e| e.start < cn.unwrap().start) {
                let e = exon.unwrap();
                self.exon = Some(e);
                if let Some(c) = self.copy_number {
                    self.add_overlap(e, c);
                }
                ei += 1;
            } else {
                let c = cn.unwrap();
                self.copy_number = Some(c);
                if let Some(e) = self.exon {
                    self.add_overlap(e, c);
                }
                ci += 1;
            }
        }
        let relative = if ploidy > 0.0 { self.min_cn / ploidy } else { 0.0 };
        GeneCopyNumber {
            chromosome: gene.chromosome.clone(),
            start: gene.gene_start,
            end: gene.gene_end,
            gene_name: gene.gene_name.clone(),
            trans_name: trans.trans_name.clone(),
            is_canonical: trans.is_canonical,
            chromosome_band: gene.karyotype_band.clone(),
            max_copy_number: self.max_cn,
            min_copy_number: self.min_cn,
            min_minor_allele_copy_number: self.min_minor,
            relative_min_copy_number: relative,
            somatic_regions: self.somatic_count,
            min_regions: self.min_regions,
            min_region_start: self.min_region_start,
            min_region_end: self.min_region_end,
            depth_window_count: self.depth_window_count,
            gc_content: self.min_gc,
            min_region_start_support: self.min_start_support,
            min_region_end_support: self.min_end_support,
            min_region_method: self.min_method,
            driver_type: DriverType::Unknown,
            reported_status: ReportedStatus::None,
        }
    }
}

/// GeneCopyNumberBuilder.createGeneCopyNumbers
pub fn create_gene_copy_numbers(v38: bool, ensembl: &EnsemblDataCache, copy_numbers: &[PurpleCopyNumber], ploidy: f64) -> Vec<GeneCopyNumber> {
    let mut result = Vec::new();
    for c in Chromosome::ALL {
        let chr_string = c.versioned(v38);
        let chr_copy_numbers: Vec<&PurpleCopyNumber> = copy_numbers.iter().filter(|x| x.chromosome == chr_string).collect();
        let Some(genes) = ensembl.chr_gene_data_map().get(&chr_string) else {
            panic!("NullPointerException: no Ensembl genes for chromosome {chr_string}");
        };
        for g in genes {
            let Some(transcripts) = ensembl.transcripts(&g.gene_id) else {
                panic!("NullPointerException: no Ensembl transcripts for gene {}", g.gene_id);
            };
            for t in transcripts {
                let b = Builder {
                    copy_numbers: &chr_copy_numbers,
                    min_cn: f64::MAX,
                    min_minor: f64::MAX,
                    max_cn: -f64::MAX,
                    somatic_count: 0,
                    previous: None,
                    previous_cn: -f64::MAX,
                    exon: None,
                    copy_number: None,
                    min_regions: 0,
                    min_region_start: 0,
                    min_region_end: 0,
                    min_start_support: SegmentSupport::None,
                    min_end_support: SegmentSupport::None,
                    min_method: CopyNumberMethod::Unknown,
                    depth_window_count: 0,
                    min_gc: 0.0,
                };
                let gcn = b.create(g, t, ploidy);
                if gcn.somatic_regions > 0 {
                    result.push(gcn);
                }
            }
        }
    }
    result
}

/// GeneCopyNumberFile.write
pub fn write_gene_copy_numbers(path: &str, genes: &[GeneCopyNumber]) -> anyhow::Result<()> {
    let f = |v: f64| hmf_common::fmt::decimal_format_min(v, 4, 4);
    let mut out = String::with_capacity(genes.len() * 200 + 400);
    out.push_str("chromosome\tstart\tend\tgene\tminCopyNumber\tmaxCopyNumber\tsomaticRegions\ttranscriptId\tisCanonical\tchromosomeBand\tminRegions\tminRegionStart\tminRegionEnd\tminRegionStartSupport\tminRegionEndSupport\tminRegionMethod\tminMinorAlleleCopyNumber\trelativeMinCopyNumber\tdepthWindowCount\tgcContent\treportedStatus\tdriverType\n");
    for g in genes {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            g.chromosome,
            g.start,
            g.end,
            g.gene_name,
            f(g.min_copy_number),
            f(g.max_copy_number),
            g.somatic_regions,
            g.trans_name,
            g.is_canonical,
            g.chromosome_band,
            g.min_regions,
            g.min_region_start,
            g.min_region_end,
            g.min_region_start_support,
            g.min_region_end_support,
            g.min_region_method,
            f(g.min_minor_allele_copy_number),
            f(g.relative_min_copy_number),
            g.depth_window_count,
            f(g.gc_content),
            g.reported_status,
            g.driver_type
        ));
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}
