//! Reportability and hmf-common ReportablePredicate: the HOTSPOT and REPORTED flags (-set_reportable).

use std::collections::HashMap;

use hmf_common::driver_panel::{DriverCategory, DriverGene};

use crate::constants::{HOTSPOT_FLAG, REPORTED_FLAG};
use crate::impact::builder::{parse_alt_transcript_info, VariantImpact};
use crate::impact::effects::{effects_to_list, CodingEffect};
use crate::variant::{VariantData, VariantTier};

struct ReportablePredicate {
    genes: HashMap<String, DriverGene>,
}

impl ReportablePredicate {
    fn new(category: DriverCategory, driver_genes: &[DriverGene]) -> ReportablePredicate {
        let genes = driver_genes
            .iter()
            .filter(|g| g.likelihood_type == category && g.report_somatic())
            .map(|g| (g.gene.clone(), g.clone()))
            .collect();
        ReportablePredicate { genes }
    }

    fn is_reportable(&self, vi: &VariantImpact, is_hotspot: bool) -> bool {
        if self.is_reportable_effects(&vi.gene_name, is_hotspot, &vi.canonical_effect) {
            return true;
        }
        if vi.other_reportable_effects.is_empty() {
            return false;
        }
        parse_alt_transcript_info(&vi.other_reportable_effects).iter().any(|x| self.is_reportable_effects(&vi.gene_name, is_hotspot, &x.effects))
    }

    fn is_reportable_effects(&self, gene: &str, is_hotspot: bool, effects: &str) -> bool {
        let Some(g) = self.genes.get(gene) else { return false };
        if is_hotspot && g.report_somatic_hotspot {
            return true;
        }
        // any reportable effect present, not only the top one
        let effects: Vec<CodingEffect> = effects_to_list(effects).into_iter().map(CodingEffect::effect).collect();
        if g.report_missense_and_inframe && effects.contains(&CodingEffect::Missense) {
            return true;
        }
        if g.report_splice && effects.contains(&CodingEffect::Splice) {
            return true;
        }
        if g.report_nonsense_and_frameshift && effects.contains(&CodingEffect::NonsenseOrFrameshift) {
            return true;
        }
        false
    }
}

pub struct Reportability {
    onco: ReportablePredicate,
    tsg: ReportablePredicate,
}

impl Reportability {
    pub fn new(driver_genes: &[DriverGene]) -> Reportability {
        Reportability { onco: ReportablePredicate::new(DriverCategory::Onco, driver_genes), tsg: ReportablePredicate::new(DriverCategory::Tsg, driver_genes) }
    }

    pub fn is_reported(&self, vi: Option<&VariantImpact>, is_hotspot: bool) -> bool {
        let Some(vi) = vi else { return false };
        self.onco.is_reportable(vi, is_hotspot) || self.tsg.is_reportable(vi, is_hotspot)
    }

    pub fn set_reportability(&self, variant: &mut VariantData, vi: Option<&VariantImpact>) -> anyhow::Result<()> {
        let is_hotspot = variant.tier() == VariantTier::Hotspot;
        if is_hotspot {
            super::put_once(variant.context_mut(), HOTSPOT_FLAG, hmf_common::vcf::InfoValue::Flag)?;
        }
        if variant.filters.is_empty() && self.is_reported(vi, is_hotspot) {
            super::put_once(variant.context_mut(), REPORTED_FLAG, hmf_common::vcf::InfoValue::Flag)?;
        }
        Ok(())
    }
}
