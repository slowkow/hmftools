//! SpecificSpliceRegion and SpliceClassifier.isWithinSpecialSpliceRegion (upstream 4adc014589): splice regions
//! beyond the standard ones, for a named transcript and exon, read from the resource
//! `splice/special_splice_regions.tsv`. Java carries the file inside its jar and reads it with
//! getResourceAsStream; the port carries the same file inside its binary (`resources/splice/`).
//!
//! The only region upstream lists is MET's exon 14 (ENST00000397752, positive strand): intronic bases -30 to -10
//! before the exon's start, the 3' splice site region whose variants make exon 14 skipping. A variant there is
//! marked as in a splice region (VariantTransImpact.markSpliceRegion, and so the IMPACT field's canonical splice
//! region and the transcript file's SpliceRegion column); its effects are not changed.

use std::sync::OnceLock;

use anyhow::{anyhow, bail, Context, Result};

use hmf_common::ensembl::{ExonData, TranscriptData};
use hmf_common::regions::positions_overlap;

use crate::variant::VariantData;

/// The resource, as the jar has it
const SPECIAL_SPLICE_REGIONS_FILE: &str = include_str!("../../resources/splice/special_splice_regions.tsv");

#[derive(Debug, Clone, PartialEq)]
pub struct SpecificSpliceRegion {
    pub gene_name: String,
    pub transcript_name: String,
    pub exon_rank: i32,
    /// the exon's start (its lower position) is the boundary, else its end
    pub use_start: bool,
    /// in bases relative to the exon boundary: -30 to -10 are the 21 bases from 30 before the boundary to 10
    /// before it
    pub base_start: i32,
    pub base_end: i32,
}

/// The regions of the resource, read once (Java reads them for each ImpactClassifier)
pub fn special_splice_regions() -> &'static [SpecificSpliceRegion] {
    static REGIONS: OnceLock<Vec<SpecificSpliceRegion>> = OnceLock::new();
    REGIONS.get_or_init(|| load(SPECIAL_SPLICE_REGIONS_FILE).expect("the built-in special splice regions resource"))
}

/// SpecificSpliceRegion.load: a header with the columns GeneName, TransName, ExonRank, UseStart, BaseStart and
/// BaseEnd, then one region per line
pub fn load(text: &str) -> Result<Vec<SpecificSpliceRegion>> {
    let mut lines = text.lines();
    let header = lines.next().ok_or_else(|| anyhow!("special splice regions: no header"))?;
    let columns: Vec<&str> = header.split('\t').collect();
    let index = |name: &str| columns.iter().position(|c| *c == name).ok_or_else(|| anyhow!("special splice regions: no column {name}"));
    let (gene, trans, exon, use_start, base_start, base_end) =
        (index("GeneName")?, index("TransName")?, index("ExonRank")?, index("UseStart")?, index("BaseStart")?, index("BaseEnd")?);

    let mut regions = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let value = |i: usize| values.get(i).copied().ok_or_else(|| anyhow!("special splice regions: short line '{line}'"));
        let int = |i: usize| -> Result<i32> { value(i)?.parse::<i32>().with_context(|| format!("special splice regions: line '{line}'")) };
        if values.len() < columns.len() {
            bail!("special splice regions: short line '{line}'");
        }
        regions.push(SpecificSpliceRegion {
            gene_name: value(gene)?.to_string(),
            transcript_name: value(trans)?.to_string(),
            exon_rank: int(exon)?,
            // Boolean.parseBoolean
            use_start: value(use_start)?.eq_ignore_ascii_case("true"),
            base_start: int(base_start)?,
            base_end: int(base_end)?,
        });
    }
    Ok(regions)
}

/// SpliceClassifier.isWithinSpecialSpliceRegion: the variant's positions (Position to EndPosition: for a DEL the
/// base before the deleted ones to the base after them, for an INS the bases on either side) overlap a region
/// listed for this transcript and exon
pub fn is_within_special_splice_region(variant: &VariantData, trans: &TranscriptData, exon: &ExonData, regions: &[SpecificSpliceRegion]) -> bool {
    for region in regions {
        if region.transcript_name != trans.trans_name || region.exon_rank != exon.rank {
            continue;
        }
        let exon_boundary = if region.use_start { exon.start } else { exon.end };
        let range_start = exon_boundary + region.base_start;
        let range_end = exon_boundary + region.base_end;
        if positions_overlap(variant.position, variant.end_position, range_start, range_end) {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_resource_lists_met_exon_14() {
        let regions = special_splice_regions();
        assert_eq!(
            regions,
            &[SpecificSpliceRegion {
                gene_name: "MET".to_string(),
                transcript_name: "ENST00000397752".to_string(),
                exon_rank: 14,
                use_start: true,
                base_start: -30,
                base_end: -10,
            }]
        );
    }

    #[test]
    fn load_reads_columns_by_name() {
        let text = "UseStart\tGeneName\tBaseEnd\tTransName\tBaseStart\tExonRank\nFALSE\tGENE_A\t12\tTRANS_A\t3\t2\n";
        let regions = load(text).unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!((regions[0].use_start, regions[0].base_start, regions[0].base_end, regions[0].exon_rank), (false, 3, 12, 2));
        assert!(load("GeneName\tTransName\n").is_err());
        assert!(load("GeneName\tTransName\tExonRank\tUseStart\tBaseStart\tBaseEnd\nA\tB\tx\ttrue\t1\t2").is_err());
    }
}
