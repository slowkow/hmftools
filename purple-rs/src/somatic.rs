//! SomaticVariantCache (the PAVE somatic VCF with hotspot tagging) and the somatic enrichment written to
//! purple.somatic.vcf.gz.

use anyhow::Result;
use std::collections::BTreeMap;

use hmf_common::variant_context::{Value, VcfHeader, VcfReader, VcfRecord};

use crate::config::PurpleConfig;
use crate::reference::Hotspot;
use crate::types::*;
use crate::variant::Variant;

/// HotspotEnrichment.processVariant: HOTSPOT for an exact match, NEAR_HOTSPOT within 5 bases
pub fn hotspot_enrich(rec: &mut VcfRecord, hotspots: &BTreeMap<Chromosome, Vec<Hotspot>>) {
    if !Chromosome::contains(&rec.chrom) {
        return;
    }
    if rec.info.has("HOTSPOT") || rec.info.has("NEAR_HOTSPOT") {
        return;
    }
    let Some(list) = hotspots.get(&chr(&rec.chrom)) else { return };
    let start = rec.pos;
    let end = start + rec.ref_allele.len() as i32 - 1 + 5;
    let mut near = false;
    for h in list {
        if start == h.position && h.ref_bases == rec.ref_allele && rec.alts.iter().any(|a| *a == h.alt) {
            rec.info.put("HOTSPOT", Value::Bool(true));
            return;
        }
        let h_end = h.position + h.ref_bases.len() as i32 - 1 + 5;
        if positions_overlap(start, end, h.position, h_end) {
            near = true;
        }
        if h.position > end {
            break;
        }
    }
    if near {
        rec.info.put("NEAR_HOTSPOT", Value::Bool(true));
    }
}

/// HlaCommon.containsPosition (tumor-only mode): the HLA genes of the Ensembl cache
pub struct HlaGenes {
    pub genes: Vec<(String, i32, i32)>,
}

impl HlaGenes {
    pub fn contains(&self, chromosome: &str, position: i32) -> bool {
        if chromosome != "6" && chromosome != "chr6" {
            return false;
        }
        self.genes.iter().any(|(c, s, e)| c == chromosome && position >= *s && position <= *e)
    }
}

pub struct SomaticVariantCache {
    pub header: Option<VcfHeader>,
    pub variants: Vec<Variant>,
    pub tinc_level: f64,
    pub reference_ordinal: i32,
    pub tumor_ordinal: i32,
    pub snp_count: i32,
    pub indel_count: i32,
}

impl SomaticVariantCache {
    pub fn empty() -> SomaticVariantCache {
        SomaticVariantCache { header: None, variants: Vec::new(), tinc_level: 0.0, reference_ordinal: -1, tumor_ordinal: -1, snp_count: 0, indel_count: 0 }
    }

    pub fn has_data(&self) -> bool {
        self.header.is_some()
    }

    /// SomaticVariantCache.loadSomatics
    pub fn load(&mut self, path: &str, config: &PurpleConfig, hotspots: &BTreeMap<Chromosome, Vec<Hotspot>>, hla: Option<&HlaGenes>) -> Result<()> {
        if path.is_empty() {
            return Ok(());
        }
        let mut reader = VcfReader::open(path)?;
        let header = reader.header.clone();
        // SageVcfTags.parseTincLevel
        self.tinc_level = header.meta_value("tincLevel").and_then(|v| hmf_common::java_text::parse_double(&v)).unwrap_or(0.0);
        let o = crate::sv_cache::genotype_ordinals(&header.samples, config.reference_id.as_deref(), Some(&config.tumor_id));
        self.reference_ordinal = o.reference;
        self.tumor_ordinal = o.tumor;
        let tumor_only = config.tumor_only_mode();
        while let Some(rec) = reader.next_record()? {
            let mut v = Variant::new(rec, &header, Some(&config.tumor_id), config.reference_id.as_deref());
            if tumor_only && hla.is_some_and(|h| h.contains(v.chromosome(), v.position())) {
                continue;
            }
            if config.exclude_on_specific_region(v.chromosome(), v.position()) {
                continue;
            }
            hotspot_enrich(&mut v.rec, hotspots);
            if v.is_pass() {
                if v.vtype == VariantType::Indel {
                    self.indel_count += 1;
                } else {
                    self.snp_count += 1;
                }
            }
            self.variants.push(v);
        }
        self.header = Some(header);
        crate::info!("loaded {} somatic variants from {}", self.variants.len(), path);
        Ok(())
    }
}
