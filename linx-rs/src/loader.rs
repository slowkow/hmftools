//! SvFileLoader: the SVs of PURPLE's somatic or germline SV VCF.

use anyhow::Result;
use hmf_common::sv::{
    always_pass, read_vcf, StructuralVariant, StructuralVariantData, StructuralVariantFactory, VcfRecord, ALLELE_FRACTION, INFERRED,
    PASS_FILTER, PON_COUNT, PON_FILTER_PON,
};

/// loadSomaticVariantsFromVcf: PURPLE's enriched SVs, every record kept
pub fn load_somatic_variants(vcf_file: &str) -> Result<Vec<StructuralVariantData>> {
    let vcf = read_vcf(vcf_file)?;
    let mut factory = StructuralVariantFactory::build(always_pass);
    factory.set_enrich(true);
    for r in vcf.records {
        factory.add_variant_context(r)?;
    }
    let out: Vec<StructuralVariantData> =
        factory.results().iter().enumerate().map(|(i, v)| StructuralVariantData::convert(v, i as i32)).collect();
    hmf_common::info!("loaded {} SV data records from VCF file: {}", out.len(), vcf_file);
    Ok(out)
}

/// GermlineFilter: no filter, PASS or PON alone
fn germline_filter(r: &VcfRecord) -> bool {
    if r.filters.is_empty() {
        return true;
    }
    if r.filters.len() > 1 {
        return false;
    }
    r.filters.iter().any(|x| x == PASS_FILTER || x == PON_FILTER_PON)
}

/// loadGermlineVariantsFromVcf
pub fn load_germline_variants(vcf_file: &str, sample_id: &str) -> Result<Vec<StructuralVariantData>> {
    let vcf = read_vcf(vcf_file)?;
    let Some(tumor_ordinal) = vcf.samples.iter().position(|s| s == sample_id) else {
        hmf_common::error!("vcf({}) missing sampleId(sampleId)", vcf_file);
        return Ok(Vec::new());
    };
    let reference_ordinal = if vcf.samples.len() > 1 { if tumor_ordinal == 0 { 1 } else { 0 } } else { -1 };
    let mut factory = StructuralVariantFactory::build(germline_filter);
    factory.set_genotype_ordinals(reference_ordinal, tumor_ordinal as i32);
    for r in vcf.records {
        factory.add_variant_context(r)?;
    }
    let mut out = Vec::new();
    for (i, v) in factory.results().iter().enumerate() {
        out.push(convert_germline_sv_data(v, i as i32)?);
    }
    hmf_common::info!("loaded {} germline SV data records from VCF file: {}", out.len(), vcf_file);
    Ok(out)
}

/// convertGermlineSvData: an assumed copy number of 2, the allele fraction of the first genotype
fn convert_germline_sv_data(var: &StructuralVariant, sv_id: i32) -> Result<StructuralVariantData> {
    const ASSUMED_CN: f64 = 2.0;
    let mut d = StructuralVariantData::convert(var, sv_id);
    let af_start = var.start_context.genotype_attribute_as_double(0, ALLELE_FRACTION, 0.0)?;
    d.junction_copy_number = ASSUMED_CN;
    d.start_af = af_start;
    d.adjusted_start_af = af_start;
    d.adjusted_start_copy_number = ASSUMED_CN;
    d.adjusted_start_copy_number_change = ASSUMED_CN - af_start;
    d.pon_count = var.start_context.attribute_as_int(PON_COUNT, 0)?;
    if let Some(end) = &var.end_context {
        let af_end = end.genotype_attribute_as_double(0, ALLELE_FRACTION, 0.0)?;
        d.end_af = af_end;
        d.adjusted_end_af = af_end;
        d.adjusted_end_copy_number = ASSUMED_CN;
        d.adjusted_end_copy_number_change = ASSUMED_CN - af_end;
    } else {
        d.end_af = 0.0;
        d.adjusted_end_af = 0.0;
        d.adjusted_end_copy_number = 0.0;
        d.adjusted_end_copy_number_change = 0.0;
    }
    Ok(d)
}

/// createSvData's filter: unfiltered, PASS or INFERRED (every SV in germline mode)
pub fn keep_sv(filter: &str, is_germline: bool) -> bool {
    filter.is_empty() || filter == PASS_FILTER || filter == INFERRED || is_germline
}
