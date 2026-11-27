//! GermlineSvCache: the germline SV VCF, the TINC transfer of germline SVs to the somatic set, the
//! copy-number annotation of the germline SVs and purple.sv.germline.vcf.gz.

use anyhow::Result;

use hmf_common::variant_context::{EncodeOptions, VcfHeader, VcfReader, VcfRecord, VcfWriter};

use crate::constants::*;
use crate::copynumber::PurpleCopyNumber;
use crate::region::ObservedRegion;
use crate::sv::{GenotypeOrdinals, StructuralVariant, SvCollection};
use crate::sv_cache::{add_enriched_contexts, genotype_ordinals, sv_output_header, SomaticSvCache};
use crate::types::*;

const SV_GERMLINE_AF_THRESHOLD: f64 = 0.1;
const SV_GERMLINE_AD_THRESHOLD: f64 = 0.01;

pub struct GermlineSvCache {
    collection: Option<SvCollection>,
    ordinals: GenotypeOrdinals,
    variants: Vec<StructuralVariant>,
    modified: bool,
}

/// getOrCalculateAlleleFrequency
fn allele_frequency(r: &VcfRecord, ordinal: i32) -> f64 {
    let g = r.decoded_genotype_at(ordinal as usize).expect("genotype ordinal (decoded by the SV factory)");
    if g.attrs.has("AF") {
        return g.attr_double("AF", 0.0);
    }
    let total = g.attr_int("REF", 0) + g.attr_int("REFPAIR", 0);
    let frags = g.attr_int("VF", 0);
    frags as f64 / (frags + total) as f64
}

impl GermlineSvCache {
    pub fn empty() -> GermlineSvCache {
        GermlineSvCache { collection: None, ordinals: GenotypeOrdinals { reference: -1, tumor: -1 }, variants: Vec::new(), modified: false }
    }

    pub fn load(version: &str, input: &str, reference_id: Option<&str>, tumor_id: &str) -> Result<GermlineSvCache> {
        if input.is_empty() {
            return Ok(Self::empty());
        }
        let mut reader = VcfReader::open(input)?;
        let header = sv_output_header(version, &reader.header);
        let ordinals = genotype_ordinals(&header.samples, reference_id, Some(tumor_id));
        let mut collection = SvCollection::new(header);
        // the genotypes are decoded where the SV factory reads them (VariantContextCollection.variants)
        while let Some(r) = reader.next_record()? {
            collection.add(r);
        }
        let mut cache = GermlineSvCache { collection: Some(collection), ordinals, variants: Vec::new(), modified: true };
        let n = cache.germline_variants().len();
        crate::info!("loaded {} germline SVs from {}", n, input);
        Ok(cache)
    }

    pub fn header(&self) -> Option<&VcfHeader> {
        self.collection.as_ref().map(|c| &c.header)
    }

    /// VariantContextCollection.variants
    pub fn germline_variants(&mut self) -> &[StructuralVariant] {
        if self.modified {
            self.modified = false;
            self.variants = match self.collection.as_mut() {
                Some(c) => c.variants(),
                None => Vec::new(),
            };
        }
        &self.variants
    }

    fn records(&self) -> &[VcfRecord] {
        self.collection.as_ref().map(|c| c.records.as_slice()).unwrap_or(&[])
    }

    /// checkTincVariants: unfiltered germline SVs that pass the germline filters once the tumor
    /// contamination of the normal is allowed for move to the somatic set
    pub fn check_tinc_variants(&mut self, somatic: &mut SomaticSvCache, tinc_level: f64) {
        if tinc_level == 0.0 || self.collection.is_none() {
            return;
        }
        let variants = self.germline_variants().to_vec();
        let tinc_adjustment = GERMLINE_SV_TINC_FACTOR * tinc_level + GERMLINE_SV_TINC_MARGIN;
        let mut transferred = vec![false; variants.len()];
        let mut n_transferred = 0;
        for (vi, v) in variants.iter().enumerate() {
            let records = self.records();
            let ctx = &records[v.start_context];
            if ctx.is_filtered() || ctx.info.has("SAGA") {
                continue;
            }
            let rg = ctx.decoded_genotype_at(self.ordinals.reference as usize).expect("reference genotype");
            let ref_non_alt = rg.attr_int("REF", 0) + rg.attr_int("REFPAIR", 0);
            let ref_alt = rg.attr_int("VF", 0);
            let ref_depth = (ref_alt + ref_non_alt) as f64;
            let ref_af = ref_alt as f64 / ref_depth;
            let tg = ctx.decoded_genotype_at(self.ordinals.tumor as usize).expect("tumor genotype");
            let tumor_alt = tg.attr_int("VF", 0);
            let tumor_af = allele_frequency(ctx, self.ordinals.tumor);
            let mut recovery = tumor_af * tinc_adjustment;
            if v.hotspot {
                recovery *= GERMLINE_SV_TINC_HOTSPOT_MULTIPLIER;
            }
            let adjusted_ref_af = ref_af - recovery;
            let adjusted_ref_alt = jround(ref_alt as f64 - ref_depth * recovery) as i32 as f64;
            let ad_ratio = adjusted_ref_alt / tumor_alt.max(1) as f64;
            if adjusted_ref_af >= SV_GERMLINE_AF_THRESHOLD * tumor_af && ad_ratio >= SV_GERMLINE_AD_THRESHOLD {
                continue;
            }
            let pon_count = ctx.info.as_int("PON_COUNT", 0);
            let mut recs = vec![ctx.clone()];
            if let Some(e) = v.end_context {
                recs.push(records[e].clone());
            }
            somatic.add_tinc_records(recs, pon_count);
            transferred[vi] = true;
            n_transferred += 1;
        }
        if n_transferred > 0 {
            crate::info!("tincLevel({}) recovered {} germline SVs", hmf_common::fmt::java_double_to_string(tinc_level), n_transferred);
            let collection = self.collection.as_mut().unwrap();
            let old = std::mem::take(&mut collection.records);
            for (vi, v) in variants.iter().enumerate() {
                if transferred[vi] {
                    continue;
                }
                collection.add(old[v.start_context].clone());
                if let Some(e) = v.end_context {
                    collection.add(old[e].clone());
                }
            }
            self.modified = true;
        }
    }

    /// annotateCopyNumberInfo: the germline SVs rebuilt with PURPLE_AF, PURPLE_CN, PURPLE_CN_CHANGE and
    /// PURPLE_JCN from the fitted regions and copy numbers
    pub fn annotate_copy_number_info(&mut self, fitted: &[ObservedRegion], copy_numbers: &[PurpleCopyNumber], purity: f64) {
        if self.collection.is_none() {
            return;
        }
        let variants = self.germline_variants().to_vec();
        let ordinals = self.ordinals;
        let collection = self.collection.as_mut().unwrap();
        let old = std::mem::take(&mut collection.records);
        for v in &variants {
            let regions = match_fitted_regions(v, fitted);
            let cns = match_copy_numbers(v, copy_numbers);
            let mut legs: [(Option<f64>, Option<f64>, Option<f64>); 2] = [(None, None, None); 2];
            let mut change_set = [false, false];
            let mut jcn_total = 0.0;
            let mut leg_count = 0;
            for se in 0..2 {
                let Some(_leg) = (if se == 0 { Some(&v.start) } else { v.end.as_ref() }) else { continue };
                let ctx = if se == 0 { &old[v.start_context] } else { &old[v.end_context.expect("end context")] };
                let (mut af, mut acn, mut change) = (0.0, 0.0, 0.0);
                if let Some(cn) = cns[se] {
                    let mut adjusted_cn = cn.average_tumor_copy_number;
                    if let Some((ri, cn_change)) = regions[se] {
                        let c = cn_change.abs();
                        change = c;
                        change_set[se] = true;
                        if fitted[ri].germline_status == GermlineStatus::Amplification {
                            adjusted_cn += c;
                        }
                    }
                    acn = adjusted_cn;
                    if ordinals.reference >= 0 && ordinals.tumor >= 0 && purity > 0.0 {
                        let ref_af = allele_frequency(ctx, ordinals.reference);
                        let tumor_af = allele_frequency(ctx, ordinals.tumor);
                        let ref_purity = 2.0 * (1.0 - purity);
                        af = (tumor_af * (ref_purity + adjusted_cn * purity) - ref_af * ref_purity) / (adjusted_cn * purity);
                        leg_count += 1;
                        jcn_total += af * adjusted_cn;
                    }
                }
                legs[se] = (Some(af), Some(acn), Some(change));
            }
            let jcn = if leg_count > 0 { jcn_total / leg_count as f64 } else { 0.0 };
            for se in 0..2 {
                if !change_set[se] && legs[se].2.is_some() {
                    legs[se].2 = Some(jcn);
                }
            }
            let end_rec = v.end_context.map(|i| &old[i]);
            let end = v.end.as_ref().map(|_| legs[1]);
            add_enriched_contexts(collection, Some(&old[v.start_context]), end_rec, legs[0], end, Some(jcn));
        }
        self.modified = true;
    }

    /// write: the unfiltered records
    pub fn write(&self, path: &str) -> Result<()> {
        let Some(c) = self.collection.as_ref() else { return Ok(()) };
        if path.is_empty() {
            return Ok(());
        }
        let mut writer = VcfWriter::create(path, &c.header, EncodeOptions { allow_missing_fields: true })?;
        for r in &c.records {
            if r.is_filtered() {
                continue;
            }
            writer.write(r)?;
        }
        writer.close()
    }
}

fn is_short_local(v: &StructuralVariant) -> bool {
    if matches!(v.sv_type, SvType::Del | SvType::Dup | SvType::Inv | SvType::Ins) {
        let length = v.position(false).unwrap() - v.position(true).unwrap();
        return length as f64 <= WINDOW_SIZE as f64 * 1.5;
    }
    false
}

/// matchFittedRegions: per leg, the index of the germline region and the reference-normalised copy
/// number change into it
fn match_fitted_regions(v: &StructuralVariant, fitted: &[ObservedRegion]) -> [Option<(usize, f64)>; 2] {
    let mut matched = [None, None];
    if is_short_local(v) {
        return matched;
    }
    let mut ri = 0usize;
    let mut chr_matched = false;
    for se in 0..2 {
        let Some(leg) = (if se == 0 { Some(&v.start) } else { v.end.as_ref() }) else { continue };
        if se == 1 && Some(leg.chromosome.as_str()) != v.chromosome(true) {
            chr_matched = false;
        }
        while ri < fitted.len() {
            let region = &fitted[ri];
            if !chr_matched && region.chromosome == leg.chromosome {
                chr_matched = true;
            }
            if chr_matched && region.start > leg.position + WINDOW_SIZE {
                break;
            }
            let idx = ri;
            ri += 1;
            if idx == 0 || region.germline_status == GermlineStatus::Unknown {
                continue;
            }
            let within = leg.chromosome == region.chromosome && position_within(leg.position, region.min_start - WINDOW_SIZE, region.max_start + WINDOW_SIZE);
            if !within {
                continue;
            }
            let mut prev = idx as i64 - 1;
            let mut prev_region = Some(&fitted[prev as usize]);
            while prev_region.unwrap().germline_status == GermlineStatus::Unknown {
                prev -= 1;
                if prev < 0 {
                    prev_region = None;
                    break;
                }
                prev_region = Some(&fitted[prev as usize]);
            }
            let Some(p) = prev_region else { continue };
            let change = region.ref_normalised_copy_number - p.ref_normalised_copy_number;
            let implied: i8 = if change > 0.0 { -1 } else { 1 };
            if implied != leg.orientation {
                continue;
            }
            let valid = matches!(region.germline_status, GermlineStatus::HomDeletion | GermlineStatus::HetDeletion | GermlineStatus::Amplification);
            matched[se] = Some((if valid { idx } else { prev as usize }, change));
            // the Java loop breaks before its increment
            ri = idx;
            break;
        }
    }
    matched
}

fn match_copy_numbers<'a>(v: &StructuralVariant, copy_numbers: &'a [PurpleCopyNumber]) -> [Option<&'a PurpleCopyNumber>; 2] {
    let mut matched: [Option<&PurpleCopyNumber>; 2] = [None, None];
    let mut ci = 0usize;
    let mut chr_matched = false;
    let matches = |leg: &crate::sv::SvLeg, c: &PurpleCopyNumber| c.chromosome == leg.chromosome && position_within(leg.position, c.start, c.end);
    for se in 0..2 {
        let Some(leg) = (if se == 0 { Some(&v.start) } else { v.end.as_ref() }) else { continue };
        if se == 1 && Some(leg.chromosome.as_str()) != v.chromosome(true) {
            chr_matched = false;
        }
        if let Some(s) = matched[0].filter(|s| matches(leg, s)) {
            matched[1] = Some(s);
            continue;
        }
        while ci < copy_numbers.len() {
            let c = &copy_numbers[ci];
            if !chr_matched && c.chromosome == leg.chromosome {
                chr_matched = true;
            }
            if matches(leg, c) {
                matched[se] = Some(c);
                break;
            }
            if chr_matched && c.start > leg.position {
                break;
            }
            ci += 1;
        }
    }
    matched
}
