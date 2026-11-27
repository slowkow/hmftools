//! Synthetic records and variants for the unit tests (no sample data).

use hmf_common::variant_context::{Attributes, Genotypes, VcfHeader, VcfRecord, NO_LOG10_PERROR};

use crate::variant::Variant;

/// a record without genotypes (VariantContextBuilder("Source", contig, start, end, alleles).make())
pub fn record(chromosome: &str, pos: i32, ref_allele: &str, alt: &str) -> VcfRecord {
    VcfRecord {
        chrom: chromosome.to_string(),
        pos,
        id: ".".into(),
        ref_allele: ref_allele.into(),
        alts: vec![alt.into()],
        log10_p_error: NO_LOG10_PERROR,
        filters: None,
        info: Attributes::default(),
        genotypes: Genotypes::Decoded(Vec::new()),
    }
}

pub fn variant(chromosome: &str, pos: i32, ref_allele: &str, alt: &str) -> Variant {
    let header = VcfHeader::from_lines(Vec::new(), Vec::new());
    Variant::new(record(chromosome, pos, ref_allele, alt), &header, None, None)
}

/// a driver gene with nothing reported (tests set the fields they need)
pub fn driver_gene(gene: &str, category: hmf_common::driver_panel::DriverCategory) -> hmf_common::driver_panel::DriverGene {
    use hmf_common::driver_panel::{DriverGene, GermlineReporting};
    DriverGene {
        gene: gene.to_string(),
        report_missense_and_inframe: false,
        report_nonsense_and_frameshift: false,
        report_splice: false,
        report_deletion: false,
        report_het_deletion: false,
        report_loh: false,
        het_deletion_threshold: 0.6,
        report_disruption: false,
        report_amplification: false,
        amplification_ratio: 3.0,
        report_somatic_hotspot: false,
        likelihood_type: category,
        report_germline_variant: GermlineReporting::None,
        report_germline_hotspot: GermlineReporting::None,
        report_germline_disruption: GermlineReporting::None,
        report_germline_deletion: GermlineReporting::None,
        report_germline_amplification: false,
        additional_reported_transcripts: Vec::new(),
        report_pgx: false,
        report_high_expression: false,
        report_low_expression: false,
        report_novel_splice_junction: false,
    }
}

/// a record parsed from a VCF data line with one sample, SAMPLE
pub fn parse_record(line: &str) -> VcfRecord {
    let header = VcfHeader::from_lines(Vec::new(), vec!["SAMPLE".to_string()]);
    VcfRecord::parse(line, &header).expect("test VCF line")
}

/// PurityAdjusterTestFactory / createCopyNumber: a copy number with the given value and supports
pub fn copy_number(chromosome: &str, start: i32, end: i32, cn: f64, start_support: crate::types::SegmentSupport, end_support: crate::types::SegmentSupport) -> crate::copynumber::PurpleCopyNumber {
    crate::copynumber::PurpleCopyNumber {
        chromosome: chromosome.to_string(),
        start,
        end,
        baf_count: 0,
        average_actual_baf: 0.0,
        average_observed_baf: 0.0,
        average_tumor_copy_number: cn,
        depth_window_count: 1,
        segment_start_support: start_support,
        segment_end_support: end_support,
        method: crate::types::CopyNumberMethod::Unknown,
        gc_content: 0.0,
        min_start: start,
        max_start: start,
    }
}

#[derive(Clone, Copy, PartialEq)]
pub enum Hotspot {
    None,
    Hotspot,
    Near,
}

/// MiscTestUtils.createVariant: a SNP (A>G), MNP (AA>GG) or INDEL (AAA>A) in GENE_01 with the coding
/// effect, repeat count, hotspot status and a variant copy number of 2 * vaf (biallelic from 0.75)
pub fn somatic(vtype: crate::types::VariantType, effect: &str, repeat_count: i32, hotspot: Hotspot, vaf: f64) -> Variant {
    use crate::types::VariantType;
    use hmf_common::variant_context::Value;
    let (r, a) = match vtype {
        VariantType::Snp => ("A", "G"),
        VariantType::Mnp => ("AA", "GG"),
        _ => ("AAA", "A"),
    };
    let mut rec = record("1", 100, r, a);
    rec.info.put("REP_C", Value::Int(repeat_count));
    match hotspot {
        Hotspot::Hotspot => rec.info.put("HOTSPOT", Value::Bool(true)),
        Hotspot::Near => rec.info.put("NEAR_HOTSPOT", Value::Bool(true)),
        Hotspot::None => {}
    }
    if 2.0 * vaf >= 1.5 {
        rec.info.put("BIALLELIC", Value::Bool(true));
    }
    rec.info.put("PURPLE_VCN", Value::Double(2.0 * vaf));
    let impact = ["GENE_01", "TRANS_01", "missense_variant", effect, "HgvsCoding", "HgvsPrtotein", "false", "", effect, "1"];
    rec.info.put("IMPACT", Value::StrList(impact.iter().map(|s| s.to_string()).collect()));
    let header = VcfHeader::from_lines(Vec::new(), Vec::new());
    Variant::new(rec, &header, None, None)
}

/// OncoDriversTest.createGeneCopyNumber
pub fn gene_copy_number(gene: &str) -> crate::gene::GeneCopyNumber {
    use crate::types::*;
    crate::gene::GeneCopyNumber {
        chromosome: "1".into(),
        start: 0,
        end: 0,
        gene_name: gene.into(),
        trans_name: String::new(),
        is_canonical: true,
        chromosome_band: String::new(),
        max_copy_number: 0.0,
        min_copy_number: 0.0,
        min_minor_allele_copy_number: 0.0,
        relative_min_copy_number: 0.0,
        somatic_regions: 1,
        min_regions: 1,
        min_region_start: 0,
        min_region_end: 0,
        depth_window_count: 0,
        gc_content: 1.0,
        min_region_start_support: SegmentSupport::None,
        min_region_end_support: SegmentSupport::None,
        min_region_method: CopyNumberMethod::Unknown,
        driver_type: DriverType::Unknown,
        reported_status: ReportedStatus::None,
    }
}

/// CobaltTestUtils.female() / male(): typical median ratios (Y counted only for the male)
pub fn cobalt_chromosomes(female: bool) -> crate::inputs::CobaltChromosomes {
    use crate::inputs::MedianRatio;
    let m = |c: &str, r: f64, n: i32| MedianRatio { chromosome: c.to_string(), median_ratio: r, count: n };
    let mut ratios: Vec<MedianRatio> = (1..=22).map(|i| m(&i.to_string(), 1.0, 1)).collect();
    if female {
        ratios.push(m("X", 1.0, 1));
    } else {
        ratios.push(m("X", 0.5, 1));
        ratios.push(m("Y", 0.5, 1000));
    }
    crate::inputs::CobaltChromosomes::new(&ratios, true)
}

/// PurpleTestUtils.createObservedRegion
pub fn observed_region(chromosome: &str, start: i32, end: i32, baf: f64, tumor_ratio: f64, status: crate::types::GermlineStatus, baf_count: i32) -> crate::region::ObservedRegion {
    crate::region::ObservedRegion::new(chromosome, start, end, true, crate::types::SegmentSupport::None, baf_count, baf, 1, tumor_ratio, 1.0, 1.0, status, false, 0.93, start, start)
}

/// FittingConfig with every default (the options registered and none given)
pub fn default_fitting_config() -> crate::config::FittingConfig {
    let mut cb = hmf_common::config_builder::ConfigBuilder::new("Purple", "4.5");
    crate::config::register(&mut cb);
    cb.parse(&["-ensembl_data_dir".to_string(), std::env::temp_dir().to_string_lossy().into_owned()]).expect("default options");
    crate::config::FittingConfig::from(&cb, false)
}

/// FittingTestUtils.createDefaultFittedRegion
pub fn default_fitted_region(chromosome: &str, start: i32, end: i32) -> crate::region::ObservedRegion {
    use crate::types::*;
    let mut r = crate::region::ObservedRegion::new(chromosome, start, end, true, SegmentSupport::None, 1, 0.5, 1, 1.0, 1.0, 1.0, GermlineStatus::Diploid, false, 0.93, 0, 0);
    r.ref_normalised_copy_number = 2.0;
    r.tumor_copy_number = 2.0;
    r.tumor_baf = 0.5;
    r
}
