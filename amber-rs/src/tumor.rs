//! TumorAnalysis: the tumor reads at the germline heterozygous sites (the BAFs) and at the germline
//! homozygous sites (the contamination sites), counted in one pass over the tumor BAM.

use anyhow::Result;

use crate::config::AmberConfig;
use crate::constants::{CONTAMINATON_MIN_NORMAL_READ_DEPTH, QUAL_FILTERED_THRESHOLD};
use crate::contamination::TumorContamination;
use crate::evidence::{Counts, Evidence};
use crate::sites::{human_chromosome, ByChr};

/// TumorBAF: the tumor evidence at a germline heterozygous site with the normal's counts
#[derive(Clone, Copy, Debug)]
pub struct TumorBaf {
    pub normal_read_depth: i32,
    pub normal_ref_support: i32,
    pub normal_alt_support: i32,
    pub tumor: Evidence,
}

impl TumorBaf {
    pub fn from_normal(normal: &Evidence) -> TumorBaf {
        TumorBaf {
            normal_read_depth: normal.c.read_depth,
            normal_ref_support: normal.c.ref_support,
            normal_alt_support: normal.c.alt_support,
            tumor: normal.copy_site(),
        }
    }
    pub fn to_types(&self) -> crate::types::TumorBaf {
        crate::types::TumorBaf {
            normal_read_depth: self.normal_read_depth,
            normal_ref_support: self.normal_ref_support,
            normal_alt_support: self.normal_alt_support,
            tumor: self.tumor.to_position_evidence(),
        }
    }
}

/// AmberUtils.aboveQualFilter: under the threshold proportion of the reads filtered for quality
pub fn above_qual_filter(c: &Counts, threshold: f64) -> bool {
    if c.read_depth == 0 {
        return false;
    }
    let filtered = c.base_qual_filtered + c.map_qual_filtered + c.seq_tech_filtered;
    (filtered as f64 / c.read_depth as f64) < threshold
}

pub struct TumorAnalysis {
    /// the sites without indels in the tumor, by chromosome in the germline order
    pub bafs: ByChr<TumorBaf>,
    /// in genome order (Java: the order of a HashMap; see REPORT.md)
    pub contamination: Vec<TumorContamination>,
}

impl TumorAnalysis {
    pub fn run(config: &AmberConfig, het: &ByChr<Evidence>, hom: &ByChr<Evidence>) -> Result<TumorAnalysis> {
        crate::info!("processing tumor germline heterozygous({}) and homozygous({}) sites", het.len(), hom.len());
        let mut tumor_bafs: ByChr<TumorBaf> = het.map(TumorBaf::from_normal);
        let mut hom_tumor: ByChr<Evidence> = hom.map(|e| e.copy_site());

        // per chromosome: the het sites then the hom sites, sorted by position (stable)
        let mut order: Vec<Vec<(bool, usize)>> = vec![Vec::new(); 25];
        let mut positions: ByChr<Evidence> = ByChr::new();
        for ci in 0..25 {
            let mut entries: Vec<(i32, bool, usize)> = Vec::with_capacity(tumor_bafs.lists[ci].len() + hom_tumor.lists[ci].len());
            entries.extend(tumor_bafs.lists[ci].iter().enumerate().map(|(i, b)| (b.tumor.pos, false, i)));
            entries.extend(hom_tumor.lists[ci].iter().enumerate().map(|(i, e)| (e.pos, true, i)));
            entries.sort_by_key(|e| e.0);
            positions.lists[ci] =
                entries.iter().map(|&(_, is_hom, i)| if is_hom { hom_tumor.lists[ci][i] } else { tumor_bafs.lists[ci][i].tumor }).collect();
            order[ci] = entries.iter().map(|&(_, is_hom, i)| (is_hom, i)).collect();
        }
        crate::bam::process_bam(config, config.tumor_bam.as_deref().unwrap(), &mut positions)?;
        for ci in 0..25 {
            for (k, &(is_hom, i)) in order[ci].iter().enumerate() {
                let counted = positions.lists[ci][k].c;
                if is_hom {
                    hom_tumor.lists[ci][i].c = counted;
                } else {
                    tumor_bafs.lists[ci][i].tumor.c = counted;
                }
            }
        }

        let mut bafs: ByChr<TumorBaf> = ByChr::new();
        for b in tumor_bafs.values() {
            if b.tumor.c.indel_count == 0 {
                bafs.put(human_chromosome(b.tumor.chr), *b);
            }
        }
        let mut contamination = Vec::new();
        for (ci, normals) in hom.chromosomes() {
            for (normal, tumor) in normals.iter().zip(&hom_tumor.lists[ci]) {
                if normal.c.ref_support <= CONTAMINATON_MIN_NORMAL_READ_DEPTH {
                    continue;
                }
                if tumor.c.alt_support > 0 && above_qual_filter(&tumor.c, QUAL_FILTERED_THRESHOLD) {
                    contamination.push(TumorContamination {
                        chr: normal.chr,
                        pos: normal.pos,
                        ref_: normal.ref_,
                        alt: normal.alt,
                        normal: normal.c,
                        tumor: tumor.c,
                    });
                }
            }
        }
        Ok(TumorAnalysis { bafs, contamination })
    }
}
