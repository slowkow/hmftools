//! Package `visualiser.file`: the data LINX writes for the visualiser (vis_*.tsv).

use anyhow::Result;
use hmf_common::fmt::{decimal_format_fixed, java_fixed};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::write_lines;

use crate::analysis::classification::is_filtered_resolved_type;
use crate::cn::SvCnData;
use crate::fusion::ENHANCER_TRANSCRIPT_PREFIX;
use crate::gene::GeneCache;
use crate::types::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VisGeneAnnotationType {
    FUSION,
    DRIVER,
    DISRUPTION,
    PSEUDOGENE,
    EXON_LOST,
}
impl std::fmt::Display for VisGeneAnnotationType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// VisGeneData
#[derive(Clone, Debug)]
pub struct VisGeneData {
    pub cluster_id: i32,
    pub gene_id: String,
    pub gene_name: String,
    pub trans_name: String,
    pub trans_id: i32,
    pub chromosome: String,
    pub annotation_type: VisGeneAnnotationType,
    pub exon_position_offsets: JavaHashMap<i32, [i32; 2]>,
    pub exons_lost_offsets: JavaHashMap<i32, [i32; 2]>,
    /// RegionOverride: the known enhancer region (start, end) an enhancer fusion's upstream side was made from, on
    /// `chromosome` (upstream 7a4490caa6)
    pub region_override: Option<(i32, i32)>,
}

impl VisGeneData {
    pub fn new(
        cluster_id: i32, gene_id: &str, gene_name: &str, trans_name: &str, trans_id: i32, chromosome: &str, t: VisGeneAnnotationType,
    ) -> Self {
        VisGeneData {
            cluster_id,
            gene_id: gene_id.to_string(),
            gene_name: gene_name.to_string(),
            trans_name: trans_name.to_string(),
            trans_id,
            chromosome: chromosome.to_string(),
            annotation_type: t,
            exon_position_offsets: JavaHashMap::new(),
            exons_lost_offsets: JavaHashMap::new(),
            region_override: None,
        }
    }
    pub fn with_region_override(mut self, region: Option<(i32, i32)>) -> Self {
        self.region_override = region;
        self
    }
}

/// VisFusion
#[derive(Clone, Debug)]
pub struct VisFusion {
    pub cluster_id: i32,
    pub reportable: bool,
    pub gene_name_up: String,
    pub transcript_up: String,
    pub chr_up: String,
    pub pos_up: i32,
    pub strand_up: i32,
    pub region_type_up: String,
    pub fused_exon_up: i32,
    pub gene_name_down: String,
    pub transcript_down: String,
    pub chr_down: String,
    pub pos_down: i32,
    pub strand_down: i32,
    pub region_type_down: String,
    pub fused_exon_down: i32,
}

/// VisSampleData
#[derive(Default)]
pub struct VisSampleData {
    pub sample_id: String,
    pub gene_data: Vec<VisGeneData>,
    pub fusions: Vec<VisFusion>,
}

impl VisSampleData {
    pub fn set_sample_id(&mut self, s: &str) {
        self.sample_id = s.to_string();
        self.gene_data.clear();
        self.fusions.clear();
    }
    /// addGeneExonData: one entry per gene and cluster, a fusion replacing another type
    pub fn add_gene_exon_data(&mut self, g: VisGeneData) {
        for i in 0..self.gene_data.len() {
            let o = &self.gene_data[i];
            if o.gene_id == g.gene_id && o.cluster_id == g.cluster_id {
                if o.annotation_type != VisGeneAnnotationType::FUSION && g.annotation_type == VisGeneAnnotationType::FUSION {
                    self.gene_data.remove(i);
                    self.gene_data.push(g);
                }
                return;
            }
        }
        self.gene_data.push(g);
    }
}

const INFO_TYPE_NORMAL: &str = "NORMAL";
const INFO_TYPE_FOLDBACK: &str = "FOLDBACK";
const PD_FIVE_PRIME_UTR: &str = "5-Prime UTR";
const PD_THREE_PRIME_UTR: &str = "3-Prime UTR";
const PD_NON_CODING: &str = "Non Coding";

fn vis_file(dir: &str, sample: &str, kind: &str, germline: bool) -> String {
    format!("{dir}/{sample}.linx.{}{kind}.tsv", if germline { "germline." } else { "" })
}

/// VisDataWriter: the visualiser's data files for one sample
pub struct VisDataWriter {
    pub enabled: bool,
    output_dir: String,
    germline: bool,
    gene_cache: Option<&'static GeneCache>,
}

impl VisDataWriter {
    pub fn new(output_dir: &str, gene_cache: Option<&'static GeneCache>, enabled: bool, germline: bool) -> Self {
        VisDataWriter { enabled, output_dir: output_dir.to_string(), germline, gene_cache }
    }

    pub fn write_output(&self, sample: &VisSampleData, clusters: &[Cluster], variants: &[Sv], chr_cn: &JavaHashMap<String, Vec<SvCnData>>) -> Result<()> {
        if !self.enabled {
            return Ok(());
        }
        self.write_sv_data(sample, variants)?;
        self.write_segment_data(sample, clusters)?;
        self.write_gene_data(sample)?;
        self.write_copy_number_data(sample, chr_cn)?;
        self.write_fusions(sample)
    }

    fn write_sv_data(&self, sample: &VisSampleData, variants: &[Sv]) -> Result<()> {
        let mut lines = vec!["ClusterId\tChainId\tSvId\tType\tResolvedType\tIsSynthetic\tChrStart\tChrEnd\tPosStart\tPosEnd\tOrientStart\tOrientEnd\tInfoStart\tInfoEnd\tJunctionCopyNumber\tInDoubleMinute".to_string()];
        for var in variants {
            let cluster = var.cl();
            let chains = cluster.find_chains(var);
            let chain_count = if chains.is_empty() { 1 } else { chains.len() };
            let unchained_id = if chains.is_empty() { cluster.chain_id(var) } else { -1 };
            for i in 0..chain_count {
                let chain = chains.get(i);
                let chain_id = chain.map(|c| c.id()).unwrap_or(unchained_id);
                let bs = var.be(true);
                let be = var.breakend(false);
                let info = |b: &Be| if b.sv().foldback_breakend(b.uses_start()).is_none() { INFO_TYPE_NORMAL } else { INFO_TYPE_FOLDBACK };
                lines.push(
                    [
                        cluster.id().to_string(),
                        chain_id.to_string(),
                        var.id().to_string(),
                        var.sv_type().to_string(),
                        cluster.resolved_type().to_string(),
                        cluster.is_synthetic_type().to_string(),
                        bs.chromosome().to_string(),
                        be.as_ref().map(|b| b.chromosome().to_string()).unwrap_or_else(|| "-1".into()),
                        bs.position().to_string(),
                        be.as_ref().map(|b| b.position()).unwrap_or(0).to_string(),
                        bs.orientation().to_string(),
                        be.as_ref().map(|b| b.orientation()).unwrap_or(0).to_string(),
                        info(&bs).to_string(),
                        be.as_ref().map(|b| info(b).to_string()).unwrap_or_default(),
                        java_fixed(var.jcn(), 4),
                        chain.is_some_and(|c| c.is_double_minute()).to_string(),
                    ]
                    .join("\t"),
                );
            }
        }
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_sv_data", self.germline), &lines)
    }

    fn write_segment_data(&self, sample: &VisSampleData, clusters: &[Cluster]) -> Result<()> {
        let mut lines = vec!["ClusterId\tChainId\tChromosome\tPosStart\tPosEnd\tLinkPloidy\tInDoubleMinute".to_string()];
        let seg = |lines: &mut Vec<String>, cid: i32, chain_id: i32, chr: &str, s: &str, e: &str, ploidy: f64, dm: bool| {
            lines.push(format!("{cid}\t{chain_id}\t{chr}\t{s}\t{e}\t{}\t{dm}", java_fixed(ploidy, 4)));
        };
        for cluster in clusters {
            if cluster.sv_count() == 1 && is_filtered_resolved_type(cluster.resolved_type()) {
                continue;
            }
            let mut unique: Vec<Pair> = Vec::new();
            let mut chains: Vec<Chain> = Vec::new();
            for chain in cluster.chains() {
                let mut index = 0;
                while index < chains.len() {
                    if chain.jcn() > chains[index].jcn() {
                        break;
                    }
                    index += 1;
                }
                chains.insert(index, chain);
            }
            for chain in &chains {
                let starts_on_end = chain.first_sv() == chain.last_sv();
                let chain_ploidy = chain.jcn();
                if !chain.is_closed_loop() {
                    if let Some(b) = chain.open_breakend(true) {
                        seg(&mut lines, cluster.id(), chain.id(), b.chromosome(), &position_value(&b, true), &position_value(&b, false), chain_ploidy, false);
                    }
                }
                for pair in chain.linked_pairs() {
                    if unique.iter().any(|x| x.matches(&pair)) {
                        continue;
                    }
                    unique.push(pair.clone());
                    let mut link_ploidy = 0.0;
                    if cluster.requires_replication() || cluster.chains().len() > 1 {
                        for other in cluster.chains() {
                            let repeats = other.linked_pairs().iter().filter(|x| x.matches(&pair)).count() as i32;
                            link_ploidy += repeats as f64 * other.jcn();
                        }
                    } else {
                        link_ploidy = chain_ploidy;
                    }
                    let bs = pair.breakend(true);
                    let be = pair.breakend(false);
                    seg(
                        &mut lines,
                        cluster.id(),
                        chain.id(),
                        bs.chromosome(),
                        &bs.position().to_string(),
                        &be.position().to_string(),
                        link_ploidy,
                        chain.is_double_minute(),
                    );
                }
                if !chain.is_closed_loop() {
                    if let Some(b) = chain.open_breakend(false) {
                        if !starts_on_end {
                            seg(&mut lines, cluster.id(), chain.id(), b.chromosome(), &position_value(&b, true), &position_value(&b, false), chain_ploidy, false);
                        }
                    }
                }
            }
            for var in cluster.unlinked_svs() {
                let chain_id = cluster.chain_id(&var);
                for s in [SE_START, SE_END] {
                    let Some(b) = var.breakend(is_start(s)) else { continue };
                    seg(&mut lines, cluster.id(), chain_id, b.chromosome(), &position_value(&b, true), &position_value(&b, false), var.jcn(), false);
                }
            }
        }
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_segments", self.germline), &lines)
    }

    pub(crate) fn write_gene_data(&self, sample: &VisSampleData) -> Result<()> {
        let mut exons = vec!["ClusterId\tGene\tTranscript\tChromosome\tAnnotationType\tExonRank\tExonStart\tExonEnd".to_string()];
        let mut proteins = vec!["ClusterId\tTranscript\tChromosome\tStart\tEnd\tInfo".to_string()];
        let cache = self.gene_cache.expect("Ensembl data for the visualiser data");
        for g in &sample.gene_data {
            if g.cluster_id < 0 {
                continue;
            }
            // checkAddEnhancerExonRegions: an enhancer fusion's upstream side is drawn as one exon over the known enhancer
            // region it was made from (until upstream 7a4490caa6: only "@IG" transcripts, over ImmuneRegions' IG locus)
            if g.annotation_type == VisGeneAnnotationType::FUSION && g.trans_name.starts_with(ENHANCER_TRANSCRIPT_PREFIX) {
                if let Some((s, e)) = g.region_override {
                    exons.push(format!("{}\t{}\t{}\t{}\t{}\t1\t{s}\t{e}", g.cluster_id, g.trans_name, g.trans_name, g.chromosome, g.annotation_type));
                    continue;
                }
            }
            let Some(t) = cache.transcript_data(&g.gene_id, &g.trans_name) else { continue };
            if t.exons.is_empty() {
                continue;
            }
            let exon_line = |cid: i32, at: VisGeneAnnotationType, rank: i32, s: i32, e: i32| {
                format!("{cid}\t{}\t{}\t{}\t{at}\t{rank}\t{s}\t{e}", g.gene_name, t.trans_name, g.chromosome)
            };
            for exon in &t.exons {
                if g.exon_position_offsets.is_empty() {
                    exons.push(exon_line(g.cluster_id, g.annotation_type, exon.rank, exon.start, exon.end));
                } else {
                    let off = g.exon_position_offsets.get(&exon.rank).expect("exon position offsets");
                    let lost = g.exons_lost_offsets.get(&exon.rank);
                    let s = exon.start + off[SE_START];
                    let e = exon.end + off[SE_END];
                    exons.push(exon_line(g.cluster_id, g.annotation_type, exon.rank, s, e));
                    if let Some(l) = lost {
                        exons.push(exon_line(g.cluster_id, VisGeneAnnotationType::EXON_LOST, exon.rank, s + l[SE_START], e + l[SE_END]));
                    }
                }
            }
            let trans_id = if g.trans_id > 0 { g.trans_id } else { t.trans_id };
            let pd = |s: i32, e: i32, info: &str| format!("{}\t{}\t{}\t{s}\t{e}\t{info}", g.cluster_id, t.trans_name, g.chromosome);
            if let Some(list) = cache.ensembl.transcript_protein_data_map().get(&trans_id) {
                for p in list {
                    let pos = hmf_common::ensembl::get_protein_domain_positions(p, t);
                    if let (Some(s), Some(e)) = (pos[0], pos[1]) {
                        proteins.push(pd(s, e, &p.hit_description));
                    }
                }
            }
            if let (Some(cs), Some(ce)) = (t.coding_start, t.coding_end) {
                let fs = if t.strand == 1 { t.trans_start } else { ce + 1 };
                let fe = if t.strand == 1 { cs - 1 } else { t.trans_end };
                let ts = if t.strand == 1 { ce + 1 } else { t.trans_start };
                let te = if t.strand == 1 { t.trans_end } else { cs - 1 };
                if fs < fe {
                    proteins.push(pd(fs, fe, PD_FIVE_PRIME_UTR));
                }
                if ts < te {
                    proteins.push(pd(ts, te, PD_THREE_PRIME_UTR));
                }
            } else {
                proteins.push(pd(t.trans_start, t.trans_end, PD_NON_CODING));
            }
        }
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_gene_exon", self.germline), &exons)?;
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_protein_domain", self.germline), &proteins)
    }

    pub fn write_copy_number_data(&self, sample: &VisSampleData, chr_cn: &JavaHashMap<String, Vec<SvCnData>>) -> Result<()> {
        let mut lines = vec!["Chromosome\tStart\tEnd\tCopyNumber\tBAF".to_string()];
        for (chr, list) in chr_cn.iter() {
            for c in list {
                lines.push(format!(
                    "{chr}\t{}\t{}\t{}\t{}",
                    c.start_pos,
                    c.end_pos,
                    decimal_format_fixed(c.copy_number, 4),
                    decimal_format_fixed(c.actual_baf, 4)
                ));
            }
        }
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_copy_number", self.germline), &lines)
    }

    fn write_fusions(&self, sample: &VisSampleData) -> Result<()> {
        write_lines(&vis_file(&self.output_dir, &sample.sample_id, "vis_fusion", self.germline), &vis_fusion_lines(&sample.fusions))
    }

    /// the empty files of a sample without SVs
    pub fn write_empty(&self, sample_id: &str) -> Result<()> {
        let d = &self.output_dir;
        let g = self.germline;
        write_lines(&vis_file(d, sample_id, "vis_sv_data", g), &["ClusterId\tChainId\tSvId\tType\tResolvedType\tIsSynthetic\tChrStart\tChrEnd\tPosStart\tPosEnd\tOrientStart\tOrientEnd\tInfoStart\tInfoEnd\tJunctionCopyNumber\tInDoubleMinute".to_string()])?;
        write_lines(&vis_file(d, sample_id, "vis_gene_exon", g), &["ClusterId\tGene\tTranscript\tChromosome\tAnnotationType\tExonRank\tExonStart\tExonEnd".to_string()])?;
        write_lines(&vis_file(d, sample_id, "vis_segments", g), &["ClusterId\tChainId\tChromosome\tPosStart\tPosEnd\tLinkPloidy\tInDoubleMinute".to_string()])?;
        write_lines(&vis_file(d, sample_id, "vis_fusion", g), &vis_fusion_lines(&[]))?;
        write_lines(&vis_file(d, sample_id, "vis_protein_domain", g), &["ClusterId\tTranscript\tChromosome\tStart\tEnd\tInfo".to_string()])
    }
}

pub fn vis_fusion_lines(fusions: &[VisFusion]) -> Vec<String> {
    let mut lines = vec!["ClusterId\tReportable\tGeneNameUp\tTranscriptUp\tChrUp\tPosUp\tStrandUp\tRegionTypeUp\tFusedExonUp\tGeneNameDown\tTranscriptDown\tChrDown\tPosDown\tStrandDown\tRegionTypeDown\tFusedExonDown".to_string()];
    for f in fusions {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.cluster_id,
            f.reportable,
            f.gene_name_up,
            f.transcript_up,
            f.chr_up,
            f.pos_up,
            f.strand_up,
            f.region_type_up,
            f.fused_exon_up,
            f.gene_name_down,
            f.transcript_down,
            f.chr_down,
            f.pos_down,
            f.strand_down,
            f.region_type_down,
            f.fused_exon_down
        ));
    }
    lines
}

fn position_value(b: &Be, is_chain_end: bool) -> String {
    let p = b.position().to_string();
    if b.forward_orient() && b.arm() == Arm::P {
        if is_chain_end { DIRECTION_TELOMERE.to_string() } else { p }
    } else if b.reverse_orient() && b.arm() == Arm::P {
        if is_chain_end { p } else { DIRECTION_CENTROMERE.to_string() }
    } else if b.reverse_orient() && b.arm() == Arm::Q {
        if is_chain_end { p } else { DIRECTION_TELOMERE.to_string() }
    } else if is_chain_end {
        DIRECTION_CENTROMERE.to_string()
    } else {
        p
    }
}
