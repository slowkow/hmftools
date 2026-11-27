//! SomaticSvCache: the somatic SV VCF read into a sorted collection, breakends inferred at copy-number
//! changes without an SV, and the annotated purple.sv.vcf.gz (PURPLE_AF, PURPLE_CN, PURPLE_CN_CHANGE,
//! PURPLE_JCN).

use anyhow::Result;

use hmf_common::variant_context::{self as vcf, EncodeOptions, Value, VcfHeader, VcfReader, VcfRecord, VcfWriter};

use crate::config::PurpleConfig;
use crate::copynumber::{by_chromosome, enrich_svs, PurpleCopyNumber};
use crate::fitting::PurityAdjuster;
use crate::sv::{inferred_record, GenotypeOrdinals, StructuralVariant, SvCollection, SvFactory};
use crate::types::*;

/// GenotypeIds.fromVcfHeader: (reference ordinal, tumor ordinal), -1 where absent; exact names first,
/// then names containing the id
pub fn genotype_ordinals(samples: &[String], reference: Option<&str>, tumor: Option<&str>) -> GenotypeOrdinals {
    let (mut t, mut r) = (-1i32, -1i32);
    for exact in [true, false] {
        for (i, s) in samples.iter().enumerate() {
            let tm = tumor.is_some_and(|id| if exact { s == id } else { s.contains(id) });
            let rm = reference.is_some_and(|id| if exact { s == id } else { s.contains(id) });
            if t < 0 && tm {
                t = i as i32;
            } else if r < 0 && rm {
                r = i as i32;
            }
        }
    }
    GenotypeOrdinals { reference: r, tumor: t }
}

pub fn sv_output_header(version: &str, template: &VcfHeader) -> VcfHeader {
    let mut h = template.clone();
    h.add_line(vcf::meta_line("purpleVersion", version));
    h.add_line(vcf::format_line("GT", "1", "String", "Genotype"));
    h.add_line(vcf::info_line("INFERRED", "0", "Flag", "Breakend inferred from copy number transition"));
    h.add_line(vcf::info_line("TINC_RECOVERED", "0", "Flag", "Variant recovered from germline filters by TINC detection"));
    h.add_line(vcf::filter_line("INFERRED", "Breakend inferred from copy number transition"));
    h.add_line(vcf::info_line("PURPLE_AF", ".", "Float", "Purity adjusted variant allelic frequency"));
    h.add_line(vcf::info_line("PURPLE_CN", ".", "Float", "Purity adjusted copy number surrounding variant location"));
    h.add_line(vcf::info_line("PURPLE_JCN", "1", "Float", "Purity adjusted copy number of variant junction"));
    h.add_line(vcf::info_line("PURPLE_CN_CHANGE", ".", "Float", "Purity adjusted change in copy number at each breakend"));
    h
}

pub type LegValues = (Option<f64>, Option<f64>, Option<f64>);

/// SomaticSvCache.addEnrichedVariantContexts: PURPLE_AF, PURPLE_CN and PURPLE_CN_CHANGE as lists over the
/// legs (reversed for the end breakend) and PURPLE_JCN
pub fn add_enriched_contexts(out: &mut SvCollection, start_rec: Option<&VcfRecord>, end_rec: Option<&VcfRecord>, start: LegValues, end: Option<LegValues>, jcn: Option<f64>) {
    let mut af = Vec::new();
    let mut cn = Vec::new();
    let mut change = Vec::new();
    for leg in std::iter::once(start).chain(end) {
        if let Some(v) = leg.0 {
            af.push(v);
        }
        if let Some(v) = leg.1 {
            cn.push(v);
        }
        if let Some(v) = leg.2 {
            change.push(v);
        }
    }
    let build = |base: &VcfRecord, af: &[f64], cn: &[f64], change: &[f64]| {
        let mut r = base.clone();
        if !af.is_empty() {
            r.info.put("PURPLE_AF", Value::DoubleList(af.to_vec()));
        }
        if !cn.is_empty() {
            r.info.put("PURPLE_CN", Value::DoubleList(cn.to_vec()));
        }
        if !change.is_empty() {
            r.info.put("PURPLE_CN_CHANGE", Value::DoubleList(change.to_vec()));
        }
        if let Some(j) = jcn {
            r.info.put("PURPLE_JCN", Value::Double(j));
        }
        r
    };
    if let Some(s) = start_rec {
        out.add(build(s, &af, &cn, &change));
    }
    if let Some(e) = end_rec {
        af.reverse();
        cn.reverse();
        change.reverse();
        out.add(build(e, &af, &cn, &change));
    }
}

pub struct SomaticSvCache {
    pub output_vcf: String,
    pub collection: Option<SvCollection>,
    pub ordinals: GenotypeOrdinals,
    next_var_id: usize,
    /// VariantContextCollection.variants(), recomputed after additions
    variants: Vec<StructuralVariant>,
    modified: bool,
}

impl SomaticSvCache {
    pub fn empty() -> SomaticSvCache {
        SomaticSvCache { output_vcf: String::new(), collection: None, ordinals: GenotypeOrdinals { reference: -1, tumor: -1 }, next_var_id: 0, variants: Vec::new(), modified: false }
    }

    pub fn load(version: &str, input: &str, output: &str, config: &PurpleConfig) -> Result<SomaticSvCache> {
        let mut reader = VcfReader::open(input)?;
        let header = sv_output_header(version, &reader.header);
        let ordinals = genotype_ordinals(&header.samples, config.reference_id.as_deref(), Some(&config.tumor_id));
        let mut collection = SvCollection::new(header);
        let mut filtered = 0;
        while let Some(r) = reader.next_record()? {
            if config.sv_qual_filter > 0 && r.phred_scaled_qual() < config.sv_qual_filter as f64 {
                filtered += 1;
                continue;
            }
            collection.add(r);
        }
        let mut cache = SomaticSvCache { output_vcf: output.to_string(), collection: Some(collection), ordinals, next_var_id: 0, variants: Vec::new(), modified: true };
        let n = cache.somatic_variants().len();
        crate::info!("loaded {} somatic SVs from {} {}", n, input, if filtered > 0 { format!("filter({filtered})") } else { String::new() });
        Ok(cache)
    }

    /// VariantContextCollection.variants
    pub fn somatic_variants(&mut self) -> &[StructuralVariant] {
        if self.modified {
            self.modified = false;
            self.variants = match self.collection.as_mut() {
                Some(c) => c.variants(),
                None => Vec::new(),
            };
        }
        &self.variants
    }

    /// the variants without recomputing (after somatic_variants has been called since the last change)
    pub fn current_variants(&self) -> &[StructuralVariant] {
        &self.variants
    }

    /// SomaticSvCache.passingBnd: the passing SVs other than inferred and single breakends
    pub fn passing_bnd(&mut self) -> i32 {
        self.somatic_variants().iter().filter(|v| v.filter == "PASS" && v.sv_type != SvType::Inf && v.sv_type != SvType::Sgl).count() as i32
    }

    /// SomaticSvCache.inferMissingVariant
    pub fn infer_missing_variants(&mut self, copy_numbers: &[PurpleCopyNumber]) {
        let Some(collection) = self.collection.as_mut() else { return };
        for i in 1..copy_numbers.len() {
            let cn = &copy_numbers[i];
            if cn.segment_start_support == SegmentSupport::None {
                let prev = &copy_numbers[i - 1];
                let increasing = greater_than(cn.average_tumor_copy_number, prev.average_tumor_copy_number);
                let rec = inferred_record(&cn.chromosome, cn.start, cn.min_start, cn.max_start, increasing, self.next_var_id);
                self.next_var_id += 1;
                collection.add(rec);
                self.modified = true;
            }
        }
    }

    /// SomaticSvCache.addTincVariant: germline SV records added as somatic ones
    pub fn add_tinc_records(&mut self, mut records: Vec<VcfRecord>, pon_count: i32) {
        let Some(collection) = self.collection.as_mut() else { return };
        for r in records.iter_mut() {
            r.info.put("TINC_RECOVERED", Value::Bool(true));
            if pon_count > 0 {
                let filters = r.filters.get_or_insert_with(Vec::new);
                if !filters.iter().any(|f| f == "PON") {
                    filters.push("PON".into());
                }
            }
        }
        for r in records {
            collection.add(r);
        }
        self.modified = true;
    }

    /// SomaticSvCache.write (passOnly in tumor-only mode)
    pub fn write(&mut self, adjuster: Option<&PurityAdjuster>, copy_numbers: &[PurpleCopyNumber], pass_only: bool, gender: Gender) -> Result<()> {
        let Some(collection) = self.collection.as_mut() else { return Ok(()) };
        let header = collection.header.clone();
        // StructuralVariantFactory.build() over every record (no filter), with the genotype ordinals
        let mut records = collection.records.clone();
        let mut factory = SvFactory::new(&header, Some(self.ordinals));
        for i in 0..records.len() {
            factory.add(&mut records, i);
        }
        let unmatched = factory.unmatched();
        let svs = std::mem::take(&mut factory.results);
        let mut out = SvCollection::new(header.clone());
        let default_adjuster;
        let adjuster = match adjuster {
            Some(a) => a,
            None => {
                default_adjuster = PurityAdjuster::with_ratios([0.0; 24], 1.0, 1.0);
                &default_adjuster
            }
        };
        let enriched = enrich_svs(adjuster, &by_chromosome(copy_numbers), &svs);
        for e in &enriched {
            let sv = &svs[e.index];
            if gender == Gender::Female {
                let on_y = |c: Option<&str>| c.is_some_and(|c| Chromosome::from_name(c).is_some_and(is_y));
                if on_y(sv.chromosome(true)) || (sv.end.is_some() && on_y(sv.chromosome(false))) {
                    continue;
                }
            }
            let start = (e.start.adjusted_allele_frequency, e.start.adjusted_copy_number, e.start.adjusted_copy_number_change);
            let end = e.end.as_ref().map(|l| (l.adjusted_allele_frequency, l.adjusted_copy_number, l.adjusted_copy_number_change));
            let end_rec = sv.end_context.map(|i| &records[i]);
            add_enriched_contexts(&mut out, Some(&records[sv.start_context]), end_rec, start, end, e.junction_copy_number);
        }
        for i in unmatched {
            out.add(records[i].clone());
        }
        let mut writer = VcfWriter::create(&self.output_vcf, &header, EncodeOptions { allow_missing_fields: true })?;
        for r in &out.records {
            if pass_only && r.is_filtered() {
                continue;
            }
            writer.write(r)?;
        }
        writer.close()
    }
}
