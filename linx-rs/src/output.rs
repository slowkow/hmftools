//! The sample output files (LinxSvAnnotation, LinxCluster, LinxLink, LinxBreakend, LinxFusion,
//! LinxGermlineDisruption, the version file), in the layouts of hmf-common's `linx` package.

use anyhow::Result;
use hmf_common::fmt::{decimal_format_fixed, java_double_to_string, java_fixed};
use hmf_common::purple_files::write_lines;

use crate::analysis::classification::get_cluster_category;
use crate::gene::{TranscriptCodingType, TranscriptRegionType};
use crate::types::*;
use crate::utils::doubles_round;

pub const BREAKEND_ORIENTATION_UPSTREAM: &str = "Upstream";
pub const BREAKEND_ORIENTATION_DOWNSTREAM: &str = "Downstream";

pub fn sample_file(dir: &str, sample: &str, kind: &str, germline: bool) -> String {
    format!("{dir}/{sample}.linx.{}{kind}.tsv", if germline { "germline." } else { "" })
}

/// LinxBreakend
pub struct LinxBreakend {
    pub id: i32,
    pub sv_id: i32,
    pub vcf_id: String,
    pub coords: String,
    pub is_start: bool,
    pub gene: String,
    pub transcript_id: String,
    pub canonical: bool,
    pub gene_orientation: &'static str,
    pub disruptive: bool,
    pub reported_status: String,
    pub undisrupted_copy_number: f64,
    pub region_type: TranscriptRegionType,
    pub coding_type: TranscriptCodingType,
    pub biotype: String,
    pub exonic_base_phase: i32,
    pub next_splice_exon_rank: i32,
    pub next_splice_exon_phase: i32,
    pub next_splice_distance: i32,
    pub total_exon_count: i32,
    pub exon_up: i32,
    pub exon_down: i32,
}

const BREAKEND_HEADER: &str = "id\tsvId\tvcfId\tcoords\tisStart\tgene\ttranscriptId\tcanonical\tgeneOrientation\tdisruptive\treportedStatus\tundisruptedCopyNumber\tregionType\tcodingType\tbiotype\texonicBasePhase\tnextSpliceExonRank\tnextSpliceExonPhase\tnextSpliceDistance\ttotalExonCount\texonUp\texonDown";

pub fn breakend_lines(breakends: &[LinxBreakend]) -> Vec<String> {
    let mut lines = Vec::with_capacity(breakends.len() + 1);
    lines.push(BREAKEND_HEADER.to_string());
    for b in breakends {
        lines.push(
            [
                b.id.to_string(),
                b.sv_id.to_string(),
                b.vcf_id.clone(),
                b.coords.clone(),
                b.is_start.to_string(),
                b.gene.clone(),
                b.transcript_id.clone(),
                b.canonical.to_string(),
                b.gene_orientation.to_string(),
                b.disruptive.to_string(),
                b.reported_status.clone(),
                java_fixed(b.undisrupted_copy_number, 4),
                b.region_type.to_string(),
                b.coding_type.to_string(),
                b.biotype.clone(),
                b.exonic_base_phase.to_string(),
                b.next_splice_exon_rank.to_string(),
                b.next_splice_exon_phase.to_string(),
                b.next_splice_distance.to_string(),
                b.total_exon_count.to_string(),
                b.exon_up.to_string(),
                b.exon_down.to_string(),
            ]
            .join("\t"),
        );
    }
    lines
}

/// LinxFusion
pub struct LinxFusion {
    pub five_prime_breakend_id: i32,
    pub three_prime_breakend_id: i32,
    pub name: String,
    pub reported: bool,
    pub reported_type: String,
    pub reportable_reasons: String,
    pub phased: String,
    pub likelihood: String,
    pub five_prime_vcf_id: String,
    pub three_prime_vcf_id: String,
    pub five_prime_coords: String,
    pub three_prime_coords: String,
    pub chain_length: i32,
    pub chain_links: i32,
    pub chain_terminated: bool,
    pub domains_kept: String,
    pub domains_lost: String,
    pub skipped_exons_up: i32,
    pub skipped_exons_down: i32,
    pub fused_exon_up: i32,
    pub fused_exon_down: i32,
}

const FUSION_HEADER: &str = "fivePrimeBreakendId\tthreePrimeBreakendId\tname\treported\treportedType\treportableReasons\tphased\tlikelihood\tfivePrimeVcfId\tthreePrimeVcfId\tfivePrimeCoords\tthreePrimeCoords\tchainLength\tchainLinks\tchainTerminated\tdomainsKept\tdomainsLost\tskippedExonsUp\tskippedExonsDown\tfusedExonUp\tfusedExonDown";

pub fn fusion_lines(fusions: &[LinxFusion]) -> Vec<String> {
    let mut lines = vec![FUSION_HEADER.to_string()];
    for f in fusions {
        lines.push(
            [
                f.five_prime_breakend_id.to_string(),
                f.three_prime_breakend_id.to_string(),
                f.name.clone(),
                f.reported.to_string(),
                f.reported_type.clone(),
                f.reportable_reasons.clone(),
                f.phased.clone(),
                f.likelihood.clone(),
                f.five_prime_vcf_id.clone(),
                f.three_prime_vcf_id.clone(),
                f.five_prime_coords.clone(),
                f.three_prime_coords.clone(),
                f.chain_length.to_string(),
                f.chain_links.to_string(),
                f.chain_terminated.to_string(),
                f.domains_kept.clone(),
                f.domains_lost.clone(),
                f.skipped_exons_up.to_string(),
                f.skipped_exons_down.to_string(),
                f.fused_exon_up.to_string(),
                f.fused_exon_down.to_string(),
            ]
            .join("\t"),
        );
    }
    lines
}

const SVS_HEADER: &str = "vcfIdStart\tvcfIdEnd\tsvId\tcoordsStart\tcoordsEnd\ttype\tclusterId\tclusterReason\tfragileSiteStart\tfragileSiteEnd\tisFoldback\tlineTypeStart\tlineTypeEnd\tjunctionCopyNumberMin\tjunctionCopyNumberMax\tgeneStart\tgeneEnd\tlocalTopologyIdStart\tlocalTopologyIdEnd\tlocalTopologyStart\tlocalTopologyEnd\tlocalTICountStart\tlocalTICountEnd";
const CLUSTERS_HEADER: &str = "clusterId\tcategory\tsynthetic\tresolvedType\tclusterCount\tclusterDesc";
const LINKS_HEADER: &str = "clusterId\tchainId\tchainIndex\tchainCount\tlowerSvId\tupperSvId\tlowerBreakendIsStart\tupperBreakendIsStart\tchromosome\tarm\tassembled\ttraversedSVCount\tlength\tjunctionCopyNumber\tjunctionCopyNumberUncertainty\tpseudogeneInfo\tecDna";

/// generateSvDataOutput and LinxSvAnnotation.write
pub fn write_svs(path: &str, variants: &[Sv]) -> Result<()> {
    let mut lines = vec![SVS_HEADER.to_string()];
    for var in variants {
        let Some(cluster) = var.cluster() else {
            hmf_common::error!("SV({}) not assigned to a cluster", var.pos_id());
            continue;
        };
        let arm_clusters = cluster.arm_clusters.borrow();
        let Some(acs) = cluster.find_arm_cluster_id(&var.be(true)) else {
            hmf_common::error!("start breakend({}) not assigned to arm cluster", var.be(true));
            continue;
        };
        let acs = &arm_clusters[acs];
        let ace = if !var.is_sgl_breakend() { var.breakend(false).and_then(|b| cluster.find_arm_cluster_id(&b)).map(|i| &arm_clusters[i]) } else { None };
        if var.linked_svs().is_some() && var.is_sgl_breakend() {
            continue;
        }
        let d = var.sv_data();
        lines.push(
            [
                d.vcf_id_start.clone(),
                d.vcf_id_end.clone(),
                var.id().to_string(),
                var.coords_str(true),
                var.coords_str(false),
                var.sv_type().to_string(),
                cluster.id().to_string(),
                var.cluster_reason(),
                var.is_fragile_site(true).to_string(),
                var.is_fragile_site(false).to_string(),
                var.is_foldback().to_string(),
                var.line_element(true).to_str().to_string(),
                var.line_element(false).to_str().to_string(),
                decimal_format_fixed(var.jcn_min(), 4),
                decimal_format_fixed(var.jcn_max(), 4),
                var.gene_in_breakend(true, false, false),
                var.gene_in_breakend(false, false, false),
                acs.id.to_string(),
                ace.map(|a| a.id).unwrap_or(-1).to_string(),
                acs.type_str(),
                ace.map(|a| a.type_str()).unwrap_or_default(),
                acs.ti_count.to_string(),
                ace.map(|a| a.ti_count).unwrap_or(0).to_string(),
            ]
            .join("\t"),
        );
    }
    write_lines(path, &lines)
}

/// generateClusterOutput and LinxCluster.write
pub fn write_clusters(path: &str, clusters: &[Cluster]) -> Result<()> {
    let mut lines = vec![CLUSTERS_HEADER.to_string()];
    for c in clusters {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            c.id(),
            get_cluster_category(c),
            c.is_synthetic_type(),
            c.resolved_type(),
            c.sv_count(),
            c.desc()
        ));
    }
    write_lines(path, &lines)
}

/// generateLinksOutput and LinxLink.write
pub fn write_links(path: &str, clusters: &[Cluster]) -> Result<()> {
    let mut lines = vec![LINKS_HEADER.to_string()];
    for cluster in clusters {
        for chain in cluster.chains() {
            let chain_sv_count = chain.sv_count();
            let mut unique: Vec<Pair> = Vec::new();
            let links = chain.linked_pairs();
            let dm = chain.is_double_minute();
            for (ci, pair) in links.iter().enumerate() {
                if unique.iter().any(|x| x.matches(pair)) {
                    continue;
                }
                unique.push(pair.clone());
                let mut chain_index = ci.to_string();
                for (j, other) in links.iter().enumerate().skip(ci + 1) {
                    if other.matches(pair) {
                        chain_index.push(';');
                        chain_index.push_str(&j.to_string());
                    }
                }
                let bs = pair.breakend(true);
                let be = pair.breakend(false);
                let arm = if bs.arm() == be.arm() { bs.arm().to_string() } else { format!("{}_{}", bs.arm(), be.arm()) };
                lines.push(
                    [
                        cluster.id().to_string(),
                        chain.id().to_string(),
                        chain_index,
                        chain_sv_count.to_string(),
                        bs.sv().id().to_string(),
                        be.sv().id().to_string(),
                        bs.uses_start().to_string(),
                        be.uses_start().to_string(),
                        bs.chromosome().to_string(),
                        arm,
                        pair.is_assembled().to_string(),
                        pair.traversed_sv_count().to_string(),
                        pair.base_length().to_string(),
                        java_double_to_string(doubles_round(chain.jcn(), 4)),
                        java_double_to_string(doubles_round(chain.jcn_uncertainty(), 4)),
                        pair.exon_match_data(),
                        dm.to_string(),
                    ]
                    .join("\t"),
                );
            }
        }
    }
    write_lines(path, &lines)
}

pub fn write_empty_sample_files(dir: &str, sample: &str, germline: bool) -> Result<()> {
    write_lines(&sample_file(dir, sample, "svs", germline), &[SVS_HEADER.to_string()])?;
    write_lines(&sample_file(dir, sample, "clusters", germline), &[CLUSTERS_HEADER.to_string()])?;
    write_lines(&sample_file(dir, sample, "links", germline), &[LINKS_HEADER.to_string()])?;
    write_lines(&sample_file(dir, sample, "breakend", germline), &[BREAKEND_HEADER.to_string()])
}

pub fn write_empty_fusions(path: &str) -> Result<()> {
    write_lines(path, &[FUSION_HEADER.to_string()])
}

const GERMLINE_DISRUPTION_HEADER: &str = "svId\tvcfId\tchromosomeStart\tchromosomeEnd\tpositionStart\tpositionEnd\torientStart\torientEnd\ttype\tfilter\tevent\tqualScore\thomologySequenceStart\thomologySequenceEnd\tinexactHomologyStart\tinexactHomologyEnd\tjunctionCopyNumber\tadjustedAFStart\tadjustedAFEnd\tadjustedCopyNumberStart\tadjustedCopyNumberEnd\tadjustedCopyNumberChangeStart\tadjustedCopyNumberChangeEnd\tgermlineFragments\tgermlineReferenceFragmentsStart\tgermlineReferenceFragmentsEnd\ttumorFragments\ttumorReferenceFragmentsStart\ttumorReferenceFragmentsEnd\tinsertSequence\tinsertSequenceAlignments\tinsertSequenceRepeatClass\tinsertSequenceRepeatType\tgene\tclusterId\tclusterCount\tresolvedType\tlinkedByStart\tlinkedByEnd\tcohortFrequency";

pub fn germline_disruption_header() -> String {
    GERMLINE_DISRUPTION_HEADER.to_string()
}

/// a LinxGermlineDisruption row
pub fn germline_disruption_line(var: &Sv, filters: &str, gene_name: &str, cluster: &Cluster, pon_count: i32) -> String {
    let d = var.sv_data();
    let f = java_double_to_string;
    [
        var.id().to_string(),
        d.vcf_id_start.clone(),
        var.chromosome(true).to_string(),
        var.chromosome(false).to_string(),
        var.position(true).to_string(),
        var.position(false).to_string(),
        var.orientation(true).to_string(),
        var.orientation(false).to_string(),
        var.sv_type().to_string(),
        filters.to_string(),
        d.event.clone(),
        f(d.quality_score),
        d.start_homology_sequence.clone(),
        d.end_homology_sequence.clone(),
        d.inexact_homology_offset_start.to_string(),
        d.inexact_homology_offset_end.to_string(),
        f(d.junction_copy_number),
        f(d.adjusted_start_af),
        f(d.adjusted_end_af),
        f(d.adjusted_start_copy_number),
        f(d.adjusted_end_copy_number),
        f(d.adjusted_start_copy_number_change),
        f(d.adjusted_end_copy_number_change),
        d.start_normal_variant_fragment_count.to_string(),
        d.start_normal_reference_fragment_count.to_string(),
        d.end_normal_reference_fragment_count.to_string(),
        d.start_tumor_variant_fragment_count.to_string(),
        d.start_tumor_reference_fragment_count.to_string(),
        d.end_tumor_reference_fragment_count.to_string(),
        d.insert_sequence.clone(),
        d.insert_sequence_alignments.clone(),
        d.insert_sequence_repeat_class.clone(),
        d.insert_sequence_repeat_type.clone(),
        gene_name.to_string(),
        cluster.id().to_string(),
        cluster.sv_count().to_string(),
        cluster.resolved_type().to_string(),
        d.start_linked_by.clone(),
        d.end_linked_by.clone(),
        pon_count.to_string(),
    ]
    .join("\t")
}

/// VersionInfo.write: linx.version in the output directory
pub fn write_version(dir: &str) -> Result<()> {
    let build = option_env!("LINX_RS_BUILD_DATE").unwrap_or("unknown");
    std::fs::write(format!("{dir}/linx.version"), format!("version={}\nbuild.date={build}\n", crate::VERSION))?;
    Ok(())
}
