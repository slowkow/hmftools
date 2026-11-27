//! hmf-common CytoBands: the chromosome band of a position, from the resource files
//! `refgenome/cytoBands.37.tsv` and `cytoBands.38.tsv` (in `resources/`, embedded).
//!
//! As Java: the chromosome names are versioned (a "chr" prefix for GRCh38, none for GRCh37), a band's
//! start and end are both inclusive (BaseRegion.containsPosition; the file's starts are 0-based, so a
//! band's last position is also the next band's first, and the first band listed wins).

use std::collections::HashMap;

const CYTOBANDS_37: &str = include_str!("../resources/refgenome/cytoBands.37.tsv");
const CYTOBANDS_38: &str = include_str!("../resources/refgenome/cytoBands.38.tsv");

struct Band {
    start: i32,
    end: i32,
    name: String,
}

pub struct CytoBands {
    bands: HashMap<String, Vec<Band>>,
}

impl CytoBands {
    /// new CytoBands(refGenomeVersion)
    pub fn new(v37: bool) -> CytoBands {
        let text = if v37 { CYTOBANDS_37 } else { CYTOBANDS_38 };
        let mut bands: HashMap<String, Vec<Band>> = HashMap::new();
        for line in text.lines().skip(1) {
            let v: Vec<&str> = line.split('\t').collect();
            let chr = v[0];
            let chr = if v37 {
                chr.strip_prefix("chr").unwrap_or(chr).to_string()
            } else if chr.starts_with("chr") {
                chr.to_string()
            } else {
                format!("chr{chr}")
            };
            let band = Band { start: v[1].parse().expect("cytoband start"), end: v[2].parse().expect("cytoband end"), name: v[3].to_string() };
            bands.entry(chr).or_default().push(band);
        }
        CytoBands { bands }
    }

    /// getCytoBandName: None for an unknown chromosome or a position outside every band
    pub fn cyto_band_name(&self, chromosome: &str, position: i32) -> Option<&str> {
        self.bands.get(chromosome)?.iter().find(|b| position >= b.start && position <= b.end).map(|b| b.name.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bands() {
        let c38 = CytoBands::new(false);
        assert_eq!(c38.cyto_band_name("chr1", 1), Some("p36.33"));
        // the shared boundary goes to the first band
        assert_eq!(c38.cyto_band_name("chr1", 2300000), Some("p36.33"));
        assert_eq!(c38.cyto_band_name("chr1", 2300001), Some("p36.32"));
        assert_eq!(c38.cyto_band_name("1", 1), None);
        let c37 = CytoBands::new(true);
        assert_eq!(c37.cyto_band_name("1", 1), Some("p36.33"));
        assert_eq!(c37.cyto_band_name("chr1", 1), None);
    }
}
