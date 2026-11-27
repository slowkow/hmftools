//! VcfWriter (the output header and each annotated record) and TranscriptWriter (the per-transcript TSV).

use anyhow::Result;

use hmf_common::fmt::vcf_double;
use hmf_common::vcf::{filter_line, info_line, VcfHeader, VcfRecord};

use crate::constants::*;
use crate::impact::builder::VariantImpact;
use crate::impact::utils::codon_for_base;
use crate::impact::{CodingContext, ProteinContext, VariantTransImpact, SE_END, SE_START};
use crate::variant::VariantData;

/// Which annotations the run adds, for the header
pub struct HeaderOptions {
    pub write_detailed: bool,
    pub pon_enabled: bool,
    pub gnomad_enabled: bool,
    pub mappability_enabled: bool,
    pub clinvar_enabled: bool,
    pub blacklist_enabled: bool,
    pub set_reportable: bool,
}

/// VcfWriter.writeHeader: the input's header with PAVE's lines
pub fn output_header(input: &VcfHeader, o: &HeaderOptions) -> VcfHeader {
    let mut h = input.clone();
    h.add_meta_line(format!("##PaveVersion={}", crate::VERSION));
    h.add_meta_line(info_line(
        VAR_TRANS_IMPACT_ANNOTATION,
        ".",
        "String",
        "Transcript impact [Gene|GeneName|Transcript|Effects|SpliceRegion|HGVS.c|HGVS.p|RefSeqId|Exon|Codon]",
    ));
    h.add_meta_line(info_line(
        VAR_IMPACT,
        "10",
        "String",
        "Variant Impact [Gene, Transcript, CanonicalEffect, CanonicalCodingEffect, SpliceRegion, HgvsCodingImpact, HgvsProteinImpact, OtherReportableEffects, WorstCodingEffect, GenesAffected]",
    ));
    if o.write_detailed {
        h.add_meta_line(info_line(PROTEIN_CONTEXT, "1", "String", PROTEIN_CONTEXT_DESC));
    }
    if o.pon_enabled {
        h.add_meta_line(info_line(PON_COUNT, "1", "Integer", "Cohort frequency for variant"));
        h.add_meta_line(info_line(PON_MAX, "1", "Integer", "Max read depth in any sample with variant"));
        h.add_meta_line(info_line(PON_AVG_READS, "1", "Integer", "Average reads in sample with variant"));
        h.add_meta_line(info_line(PON_MULTI_STATUS, "1", "String", "Status of the variant across PONs"));
    }
    h.add_meta_line(filter_line(PON_ARTEFACT_FILTER, "Filter PON artefact"));
    h.add_meta_line(filter_line(PON_FILTER, "Filter PON variant"));
    if o.gnomad_enabled {
        h.add_meta_line(info_line(GNOMAD_FREQ, "1", "Float", GNOMAD_FREQ_DESC));
    }
    h.add_meta_line(filter_line(PON_GNOMAD_FILTER, "Filter Gnomad PON"));
    if o.mappability_enabled {
        h.add_meta_line(info_line(MAPPABILITY, "1", "Float", MAPPABILITY_DESC));
    }
    if o.clinvar_enabled {
        h.add_meta_line(info_line(CLNSIG, ".", "String", CLNSIG_DESC));
        h.add_meta_line(info_line(CLNSIGCONF, ".", "String", CLNSIGCONF_DESC));
    }
    if o.blacklist_enabled {
        h.add_meta_line(info_line(BLACKLIST_BED_FLAG, "0", "Flag", BLACKLIST_BED_FLAG));
        h.add_meta_line(info_line(BLACKLIST_VCF_FLAG, "0", "Flag", BLACKLIST_VCF_DESC));
    }
    if o.set_reportable {
        h.add_meta_line(info_line(HOTSPOT_FLAG, "0", "Flag", HOTSPOT_DESCRIPTION));
        h.add_meta_line(info_line(REPORTED_FLAG, "0", "Flag", REPORTED_DESC));
    }
    h
}

/// writeProteinInfoString: PR_CXT from a canonical impact's protein context
fn protein_info_string(pc: &ProteinContext) -> String {
    format!(
        "{},{},{},{},{},{},{},{}_{}",
        pc.ref_codon_bases_extended,
        pc.alt_codon_bases_complete,
        pc.ref_amino_acids,
        pc.alt_amino_acids,
        pc.net_ref_amino_acids,
        pc.net_alt_amino_acids,
        pc.codon_index,
        pc.net_codon_index_range[SE_START],
        pc.net_codon_index_range[SE_END]
    )
}

/// VariantTranscriptImpact.toVcfData
fn trans_impact_vcf_data(gene_name: &str, ti: &VariantTransImpact) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        ti.trans.gene_id,
        gene_name,
        ti.trans.trans_name,
        ti.effects_str(),
        ti.in_splice_region(),
        ti.hgvs_coding(),
        ti.hgvs_protein(),
        ti.trans.ref_seq_id.as_deref().unwrap_or(""),
        ti.coding.exon_rank,
        codon_for_base(ti.coding.coding_base)
    )
}

/// VcfWriter.buildVariant, in place on the variant's record: its filters with PAVE's, and PAVE's annotations
pub fn finish_record(rec: &mut VcfRecord, variant: &VariantData, vi: Option<&VariantImpact>, write_detailed: bool) {
    // the record's filters with PAVE's; PASS when there are none
    if !variant.filters.is_empty() {
        for f in &variant.filters {
            if !rec.filters.iter().any(|x| x == f) {
                rec.filters.push(f.clone());
            }
        }
        rec.filters.retain(|f| f != PASS_FILTER);
    }
    if rec.filters.is_empty() {
        rec.filters.push(PASS_FILTER.to_string());
    }

    if !variant.impacts().is_empty() {
        let mut items: Vec<String> = Vec::new();
        let mut protein_info: Option<String> = None;
        for (gene_name, impacts) in variant.impacts().iter() {
            for ti in impacts {
                items.push(trans_impact_vcf_data(gene_name, ti));
                if write_detailed && ti.trans.is_canonical && protein_info.is_none() && !rec.has_attr(PROTEIN_CONTEXT) {
                    if let Some(pc) = &ti.protein {
                        protein_info = Some(protein_info_string(pc));
                    }
                }
            }
        }
        if let Some(p) = protein_info {
            rec.put_text(PROTEIN_CONTEXT, p);
        }
        rec.put_text(VAR_TRANS_IMPACT_ANNOTATION, items.join(","));
        if let Some(vi) = vi {
            rec.put_text(VAR_IMPACT, vi.to_vcf_data().join(","));
        }
    }

    if let Some(f) = variant.gnomad_frequency() {
        if !rec.has_attr(GNOMAD_FREQ) {
            rec.put_text(GNOMAD_FREQ, vcf_double(f));
        }
    }

    if (variant.pon_sample_count > 0 || variant.pon_max_read_count > 0) && !rec.has_attr(PON_COUNT) {
        rec.put_text(PON_COUNT, variant.pon_sample_count.to_string());
        rec.put_text(PON_MAX, variant.pon_max_read_count.to_string());
    }
}

/// TranscriptWriter's header
pub fn transcript_header() -> String {
    format!(
        "{}\t{}\tGeneId\tGeneName\t{}\t{}\t{}",
        VariantData::tsv_header(),
        VariantData::extra_data_header(),
        VariantTransImpact::tsv_header(),
        CodingContext::tsv_header(),
        ProteinContext::tsv_header()
    )
}

/// TranscriptWriter.writeVariantData: a line per transcript impact of the gene
pub fn transcript_rows(variant: &VariantData, gene_name: &str, sample_id: &str, out: &mut Vec<String>) -> Result<()> {
    let Some(impacts) = variant.gene_impacts(gene_name) else { return Ok(()) };
    for ti in impacts {
        out.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            variant.tsv_data(),
            variant.extra_data_tsv(sample_id),
            ti.trans.gene_id,
            gene_name,
            ti.to_tsv(),
            ti.coding.to_tsv(),
            ti.protein.as_ref().map(|p| p.to_tsv()).unwrap_or_else(|| {
                // PA5: Java's empty protein context is ",,,0,,,0-0,,,", one field where the header has ten,
                // so these lines are short of columns. Corrected: ten tab-separated fields.
                if crate::compat::java("PA5") { ProteinContext::empty_tsv().to_string() } else { "\t\t\t0\t\t\t0-0\t\t\t".to_string() }
            })
        ));
    }
    Ok(())
}
