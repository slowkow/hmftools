//! VirusInterpreterAlgo and CoveragesAnalyzer: each virus of VIRUSBreakend's summary annotated with its
//! name, interpretation, the expected clonal coverage, and whether it is reported or blacklisted.

use anyhow::{bail, Result};

use crate::compat;
use crate::files::{AnnotatedVirus, TaxonomyDb, VirusBreakend, VirusLikelihoodType, VirusReportingDbModel};

/// CoveragesAnalyzer.calculateExpectedClonalCoverage: the tumor's mean coverage times purity over ploidy
pub fn expected_clonal_coverage(purity: f64, ploidy: f64, tumor_mean_coverage: f64) -> f64 {
    (tumor_mean_coverage * purity) / ploidy
}

/// hasAcceptablePurpleQuality: neither FAIL_NO_TUMOR nor FAIL_CONTAMINATION
pub fn has_acceptable_purple_quality(statuses: &[String]) -> bool {
    !statuses.iter().any(|s| s == "FAIL_NO_TUMOR" || s == "FAIL_CONTAMINATION")
}

pub struct VirusInterpreterAlgo<'a> {
    pub taxonomy_db: &'a TaxonomyDb,
    pub blacklisted_taxids: &'a [i32],
    pub reporting_db: &'a VirusReportingDbModel,
    pub expected_clonal_coverage: f64,
}

impl VirusInterpreterAlgo<'_> {
    pub fn analyze(&self, breakends: &[VirusBreakend], purple_qc_status: &[String]) -> Result<Vec<AnnotatedVirus>> {
        let mut out = Vec::with_capacity(breakends.len());
        for b in breakends {
            let interpretation = self.reporting_db.interpret_virus_species(b.taxid_species);
            let taxid = b.reference_taxid;
            let reported = self.report(b, self.expected_clonal_coverage, purple_qc_status);
            let blacklisted = self.blacklist(b.taxid_genus, b.taxid_species);
            if blacklisted && reported {
                // Java's message ends in the taxid after a "{}" that is never filled in (VI2)
                if compat::java("VI2") {
                    bail!("java.lang.IllegalStateException: Virus with taxid {{}} configured as reported and as blacklisted{taxid}");
                }
                bail!("virus with taxid {taxid} configured as reported and as blacklisted");
            }
            out.push(AnnotatedVirus {
                taxid,
                name: self.taxonomy_db.lookup_name(taxid),
                qc_status: b.qc_status,
                integrations: b.integrations,
                interpretation,
                percentage_covered: b.coverage,
                mean_coverage: b.mean_depth,
                expected_clonal_coverage: if has_acceptable_purple_quality(purple_qc_status) { Some(self.expected_clonal_coverage) } else { None },
                reported,
                blacklisted: Some(blacklisted),
                virus_driver_likelihood_type: self.virus_likelihood_type(b, reported),
            });
        }
        Ok(out)
    }

    pub fn virus_likelihood_type(&self, b: &VirusBreakend, reported: bool) -> VirusLikelihoodType {
        if reported { self.reporting_db.virus_likelihood_type(b.taxid_species) } else { VirusLikelihoodType::Unknown }
    }

    /// report(): a virus of the reporting database is reported when it passes VIRUSBreakend's QC and, with
    /// acceptable PURPLE quality, its covered percentage and mean depth are above the database's minimum and
    /// the expected clonal coverage (if the database gives a minimum); with unacceptable quality only if the
    /// database gives no minimum
    pub fn report(&self, b: &VirusBreakend, expected_clonal_coverage: f64, purple_qc_status: &[String]) -> bool {
        let mut reported = false;
        if self.reporting_db.has_interpretation(b.taxid_species) {
            let minimal = self.determine_minimal_coverage_virus(b.integrations, b.taxid_species);
            if has_acceptable_purple_quality(purple_qc_status) {
                match minimal {
                    Some(m) => {
                        if b.coverage > m as f64 && b.mean_depth > expected_clonal_coverage {
                            reported = true;
                        }
                    }
                    None => reported = true,
                }
            } else if minimal.is_none() {
                reported = true;
            }
        }
        reported && b.qc_status != "LOW_VIRAL_COVERAGE"
    }

    pub fn blacklist(&self, taxid_genus: i32, taxid_species: i32) -> bool {
        self.blacklisted_taxids.contains(&taxid_genus) || self.blacklisted_taxids.contains(&taxid_species)
    }

    pub fn determine_minimal_coverage_virus(&self, integrations: i32, taxid_species: i32) -> Option<i32> {
        if integrations >= 1 {
            self.reporting_db.integrated_minimal_coverage(taxid_species)
        } else {
            self.reporting_db.non_integrated_minimal_coverage(taxid_species)
        }
    }
}
