//! Blacklistings: known blacklisted regions (BED) and variants (VCF).

use anyhow::{anyhow, Context, Result};

use hmf_common::regions::ChrBaseRegion;
use hmf_common::vcf;

use crate::constants::{BLACKLIST_BED_FLAG, BLACKLIST_VCF_FLAG};
use crate::variant::VariantData;

#[derive(Default)]
pub struct Blacklistings {
    bed_regions: Vec<ChrBaseRegion>,
    vcf_entries: Vec<(String, i32, String, String)>,
}

impl Blacklistings {
    pub fn load(bed: Option<&str>, vcf_file: Option<&str>) -> Result<Blacklistings> {
        let mut b = Blacklistings::default();
        if let Some(f) = bed {
            // Files.readAllLines: every line is a region, the start 0-based
            let text = std::fs::read_to_string(f).with_context(|| format!("failed to load Blacklist BED file({f})"))?;
            for line in hmf_common::java_text::lines_of(&text) {
                let values: Vec<&str> = line.split('\t').collect();
                let bad = || anyhow!("NumberFormatException: Blacklist BED line: {line}");
                let start = values.get(1).and_then(|v| hmf_common::java_text::parse_int(v)).ok_or_else(bad)?;
                let end = values.get(2).and_then(|v| hmf_common::java_text::parse_int(v)).ok_or_else(bad)?;
                b.bed_regions.push(ChrBaseRegion::new(values[0], start + 1, end));
            }
            crate::info!("loaded {} Blacklist BED entries from file({})", b.bed_regions.len(), f);
        }
        if let Some(f) = vcf_file {
            let (header, mut reader, first) = vcf::read_header(f)?;
            let mut pending = first;
            loop {
                let line = match pending.take() {
                    Some(l) => l,
                    None => match reader.next_line()? {
                        Some(l) => String::from_utf8_lossy(l).into_owned(),
                        None => break,
                    },
                };
                if line.is_empty() {
                    continue;
                }
                let rec = vcf::VcfRecord::parse(&line, &header)?;
                if rec.alts.is_empty() {
                    anyhow::bail!("IndexOutOfBoundsException: blacklist VCF record without an alt: {line}");
                }
                b.vcf_entries.push((rec.chrom.clone(), rec.pos, rec.ref_allele.clone(), rec.alt().to_string()));
            }
            crate::info!("loaded {} blacklist VCF entries from file({})", b.vcf_entries.len(), f);
        }
        Ok(b)
    }

    pub fn enabled(&self) -> bool {
        !self.bed_regions.is_empty() || !self.vcf_entries.is_empty()
    }

    pub fn annotate(&self, variant: &mut VariantData) -> Result<()> {
        let (chr, pos) = (variant.chromosome.clone(), variant.position);
        if self.bed_regions.iter().any(|r| r.chromosome == chr && r.contains_position(pos)) {
            super::put_once(variant.context_mut(), BLACKLIST_BED_FLAG, hmf_common::vcf::InfoValue::Flag)?;
        }
        if self.vcf_entries.iter().any(|e| e.0 == chr && e.1 == pos && e.2 == variant.ref_ && e.3 == variant.alt) {
            super::put_once(variant.context_mut(), BLACKLIST_VCF_FLAG, hmf_common::vcf::InfoValue::Flag)?;
        }
        Ok(())
    }
}
