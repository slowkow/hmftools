//! CobaltScope: WholeGenome, or TargetRegions for a targeted panel (`-target_region_norm_file`): the
//! enrichment of each on-target window (TargetRegionEnrichment.loadEnrichmentFile, one TargetRegions per
//! panel label when the file has two or more labels).

use anyhow::{bail, Context, Result};
use hmf_common::java_hash::JavaHashSet;

use crate::chromosome::{Chr, PerChr};
use crate::constants::WINDOW_SIZE;

/// RefGenomeCoordinates lengths (lengths.37.tsv, lengths.38.tsv of hmf-common)
const LENGTHS_37: [i32; 24] = [
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
];
const LENGTHS_38: [i32; 24] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
];

pub enum Scope<'a> {
    WholeGenome,
    Targeted(&'a TargetRegions),
}

impl Scope<'_> {
    pub fn on_target(&self, chr: Chr, position: i32) -> Result<bool> {
        match self {
            Scope::WholeGenome => Ok(true),
            Scope::Targeted(t) => Ok(t.enrichment(chr, position)?.is_some()),
        }
    }
    /// findRegionEnrichment: 1 for the whole genome, -1 off target
    pub fn region_enrichment(&self, chr: Chr, position: i32) -> Result<f64> {
        match self {
            Scope::WholeGenome => Ok(1.0),
            Scope::Targeted(t) => Ok(t.enrichment(chr, position)?.unwrap_or(-1.0)),
        }
    }
}

/// TargetRegions: for each chromosome of the file, the enrichment of each window (None where the file has
/// no row at the window's start)
pub struct TargetRegions {
    slots: PerChr<Option<Vec<Option<f64>>>>,
}

impl TargetRegions {
    /// The enrichment of each window of each chromosome (None: off target)
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn from_slots(slots: PerChr<Option<Vec<Option<f64>>>>) -> TargetRegions {
        TargetRegions { slots }
    }

    fn enrichment(&self, chr: Chr, position: i32) -> Result<Option<f64>> {
        match self.slots.get(chr) {
            None => Ok(None),
            Some(list) => match list.get((position / WINDOW_SIZE) as usize) {
                Some(e) => Ok(*e),
                None => bail!("java.lang.IndexOutOfBoundsException: window {}:{position} beyond the reference genome length", chr.short_name()),
            },
        }
    }
}

struct Row {
    chr: Chr,
    position: i32,
    enrichment: f64,
    label: String,
}

/// TargetRegions.loadNormalisationFile
pub fn load_normalisation_file(path: &str, v38: bool) -> Result<Vec<TargetRegions>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("failed to load target-panel norm regions file({path})"))?;
    let mut lines: Vec<&str> = text.lines().collect();
    if lines.is_empty() {
        bail!("java.lang.IndexOutOfBoundsException: empty target-panel norm regions file({path})");
    }
    let header: Vec<&str> = lines.remove(0).split('\t').collect();
    let index_of = |name: &str| header.iter().position(|h| *h == name);
    let chr_index = index_of("Chromosome").unwrap_or(0);
    let pos_index = index_of("Position").unwrap_or(1);
    let enrichment_index = index_of("RelativeEnrichment").unwrap_or(2);
    let label_index = index_of("PanelLabel");
    let mut labels: JavaHashSet<String> = JavaHashSet::new();
    let mut rows = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let field = |i: usize| -> Result<&str> { values.get(i).copied().ok_or_else(|| anyhow::anyhow!("java.lang.ArrayIndexOutOfBoundsException in line: {line}")) };
        let chr_name = field(chr_index)?;
        let Some(chr) = Chr::from_name(chr_name) else { bail!("java.lang.IllegalArgumentException: No enum constant for {chr_name}") };
        let position: i32 = crate::config::java_parse_int(field(pos_index)?).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException in line: {line}"))?;
        let e = field(enrichment_index)?;
        let enrichment = if e == "NA" {
            f64::NAN
        } else {
            crate::config::java_parse_double(e).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException in line: {line}"))?
        };
        let label = match label_index {
            Some(i) => {
                let l = field(i)?.to_string();
                labels.insert(l.clone());
                l
            }
            None => String::new(),
        };
        rows.push(Row { chr, position, enrichment, label });
    }
    crate::info!("loaded {} target-panel norm regions from file({path})", rows.len());
    let lengths = if v38 { &LENGTHS_38 } else { &LENGTHS_37 };
    let build = |rows: &[&Row]| -> TargetRegions {
        let mut slots: PerChr<Option<Vec<Option<f64>>>> = PerChr::default();
        for chr in Chr::all() {
            let chr_rows: Vec<&&Row> = rows.iter().filter(|r| r.chr == chr).collect();
            if chr_rows.is_empty() {
                continue;
            }
            let count = (lengths[chr.index()] / WINDOW_SIZE) as usize;
            let mut list: Vec<Option<f64>> = vec![None; count];
            // positionToSuppliedItem: the last row at a position wins; only window starts are looked up
            for r in chr_rows {
                if r.position >= 1 && (r.position - 1) % WINDOW_SIZE == 0 {
                    let i = ((r.position - 1) / WINDOW_SIZE) as usize;
                    if i < count {
                        list[i] = Some(r.enrichment);
                    }
                }
            }
            *slots.get_mut(chr) = Some(list);
        }
        TargetRegions { slots }
    };
    if labels.len() < 2 {
        let all: Vec<&Row> = rows.iter().collect();
        return Ok(vec![build(&all)]);
    }
    crate::info!("loading normalisation for {} panels: [{}]", labels.len(), labels.iter().cloned().collect::<Vec<_>>().join(", "));
    Ok(labels
        .iter()
        .map(|label| {
            let panel: Vec<&Row> = rows.iter().filter(|r| &r.label == label).collect();
            build(&panel)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, text: &str) -> String {
        let dir = std::env::temp_dir().join(format!("cobalt_rs_targeted_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().to_string()
    }

    /// TargetRegionsTest: enrichment per window start, off target elsewhere
    #[test]
    fn enrichment_by_window() {
        let path = write("panel.tsv", "Chromosome\tPosition\tRelativeEnrichment\nchr1\t1001\t0.5\nchr1\t3001\tNA\nchr1\t3500\t2\nchr2\t1\t1.5\n");
        let t = load_normalisation_file(&path, true).unwrap();
        assert_eq!(t.len(), 1);
        let s = Scope::Targeted(&t[0]);
        let c1 = Chr::new(0);
        assert!(!s.on_target(c1, 1).unwrap());
        assert!(s.on_target(c1, 1001).unwrap());
        assert_eq!(s.region_enrichment(c1, 1001).unwrap(), 0.5);
        assert!(s.region_enrichment(c1, 3001).unwrap().is_nan());
        assert_eq!(s.region_enrichment(c1, 2001).unwrap(), -1.0);
        assert_eq!(s.region_enrichment(Chr::new(1), 1).unwrap(), 1.5);
        assert!(!s.on_target(Chr::new(2), 1).unwrap());
    }

    #[test]
    fn panels_by_label() {
        let path = write("labels.tsv", "Chromosome\tPosition\tRelativeEnrichment\tPanelLabel\nchr1\t1\t0.5\tA\nchr1\t1001\t0.7\tB\nchr1\t2001\t0.9\tA\n");
        let t = load_normalisation_file(&path, true).unwrap();
        assert_eq!(t.len(), 2);
        // HashSet<String> order: "A" (65) before "B" (66)
        assert_eq!(Scope::Targeted(&t[0]).region_enrichment(Chr::new(0), 2001).unwrap(), 0.9);
        assert_eq!(Scope::Targeted(&t[1]).region_enrichment(Chr::new(0), 1001).unwrap(), 0.7);
        assert!(!Scope::Targeted(&t[1]).on_target(Chr::new(0), 1).unwrap());
    }
}
