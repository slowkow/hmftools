//! Package `drivers`: the SV events behind PURPLE's copy number drivers (DriverGeneAnnotator,
//! AmplificationDrivers, DeletionDrivers, DriverDataCache) and LINX's homozygous disruption drivers.

use std::collections::HashMap;

use hmf_common::driver_panel::DriverGene;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::{
    driver_catalog_somatic_filename, gene_copy_number_filename, read_driver_catalog, read_gene_copy_numbers, read_purity_context,
    DriverCatalog, GeneCopyNumber,
};
use hmf_common::sv::StructuralVariantType::*;

use crate::cn::{CnDataLoader, SvCnData, TOTAL_CN_LOSS};
use crate::gene::{BreakendGene, BreakendTrans, Gene, GeneCache, Trans};
use crate::idmap::IdMap;
use crate::types::*;
use crate::utils::{copy_numbers_equal, get_chromosomal_arm, MAX_COPY_NUM_DIFF, MAX_COPY_NUM_DIFF_PERC};
use crate::vis::{VisGeneAnnotationType, VisGeneData, VisSampleData};

/// DriverEventType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DriverEventType {
    UNCLEAR,
    GAIN,
    GAIN_ARM,
    GAIN_CHR,
    LOH,
    LOH_ARM,
    LOH_CHR,
    LOH_SV_CENTRO,
    LOH_SV_TELO,
    DEL,
    HOM_DEL_DISRUPTION,
    HOM_DUP_DISRUPTION,
}
use DriverEventType as ET;

/// LinxDriver: a row of <sample>.linx.drivers.tsv
pub struct LinxDriver {
    pub cluster_id: i32,
    pub gene: String,
    pub event_type: DriverEventType,
}

pub fn linx_driver_lines(drivers: &[LinxDriver]) -> Vec<String> {
    let mut lines = vec!["clusterId\tgene\teventType".to_string()];
    for d in drivers {
        lines.push(format!("{}\t{}\t{:?}", d.cluster_id, d.gene, d.event_type));
    }
    lines
}

/// DriverGeneEvent
pub struct DriverGeneEvent {
    pub event_type: DriverEventType,
    cluster: Option<Cluster>,
    breakend_pair: [Option<Be>; 2],
    copy_number_gain: f64,
    sv_info: String,
}

impl DriverGeneEvent {
    fn new(t: DriverEventType) -> Self {
        DriverGeneEvent { event_type: t, cluster: None, breakend_pair: [None, None], copy_number_gain: 0.0, sv_info: String::new() }
    }
    pub fn cluster(&self) -> Option<Cluster> {
        if self.cluster.is_some() {
            return self.cluster.clone();
        }
        if let Some(b) = &self.breakend_pair[0] {
            return b.cluster();
        }
        self.breakend_pair[1].as_ref().and_then(|b| b.cluster())
    }
    fn add_sv_breakend_pair(&mut self, start: Option<Be>, end: Option<Be>, info: &str) {
        self.cluster = match (&start, &end) {
            (Some(s), _) => s.cluster(),
            (None, Some(e)) => e.cluster(),
            _ => self.cluster.clone(),
        };
        self.breakend_pair = [start, end];
        self.sv_info = info.to_string();
    }
}

/// GeneCopyNumberRegion
pub struct GeneCopyNumberRegion {
    pub region_start: i32,
    pub region_end: i32,
    pub min_copy_number: f64,
    pub is_germline: bool,
}

impl GeneCopyNumberRegion {
    fn from_gene_cn(g: &GeneCopyNumber) -> Self {
        GeneCopyNumberRegion { region_start: g.min_region_start, region_end: g.min_region_end, min_copy_number: g.min_copy_number, is_germline: false }
    }
}

/// calcGeneCopyNumberRegion: the exon-overlapping segment of lowest copy number
fn calc_gene_copy_number_region(trans: Trans, cn_data: Option<&Vec<SvCnData>>) -> Option<GeneCopyNumberRegion> {
    // Java dereferences the chromosome's list unchecked
    let cn_data = cn_data.expect("chromosome copy number data");
    let mut min: Option<&SvCnData> = None;
    for exon in &trans.exons {
        for c in cn_data {
            if c.end_pos < exon.start {
                continue;
            }
            if c.start_pos > exon.end {
                break;
            }
            if min.is_none_or(|m| c.copy_number < m.copy_number) {
                min = Some(c);
            }
        }
    }
    min.map(|m| GeneCopyNumberRegion { region_start: m.start_pos, region_end: m.end_pos, min_copy_number: m.copy_number, is_germline: false })
}

/// DriverGeneData
pub struct DriverGeneData {
    pub driver: DriverCatalog,
    pub gene: Gene,
    pub trans: Trans,
    pub cn_region: GeneCopyNumberRegion,
    pub arm: Arm,
    pub events: Vec<DriverGeneEvent>,
}

impl DriverGeneData {
    fn new(driver: DriverCatalog, gene: Gene, trans: Trans, cn_region: GeneCopyNumberRegion) -> Self {
        DriverGeneData { arm: get_chromosomal_arm(&gene.chromosome, gene.gene_start), driver, gene, trans, cn_region, events: Vec::new() }
    }
}

/// DriverGeneAnnotator with its DriverDataCache
pub struct DriverGeneAnnotator {
    gene_cache: &'static GeneCache,
    driver_genes: &'static JavaHashMap<String, DriverGene>,
    purple_data_path: String,
    reportable_disruption_gene_transcripts: HashMap<String, Vec<String>>,
    pub driver_catalog: Vec<DriverCatalog>,
    pub driver_gene_data: Vec<DriverGeneData>,
    gene_copy_numbers: Vec<GeneCopyNumber>,
    pub driver_output: Vec<LinxDriver>,
    sample_ploidy: f64,
    is_male: bool,
}

const MIN_AMP_PERCENT_VS_MAX: f64 = 0.33;

impl DriverGeneAnnotator {
    pub fn new(gene_cache: &'static GeneCache, driver_genes: &'static JavaHashMap<String, DriverGene>, purple_data_path: &str) -> Self {
        let mut map = HashMap::new();
        for dg in driver_genes.values() {
            if !dg.report_disruption {
                continue;
            }
            if let Some(g) = gene_cache.gene_data_by_name(&dg.gene) {
                map.insert(g.gene_id.clone(), dg.additional_reported_transcripts.clone());
            }
        }
        DriverGeneAnnotator {
            gene_cache,
            driver_genes,
            purple_data_path: purple_data_path.to_string(),
            reportable_disruption_gene_transcripts: map,
            driver_catalog: Vec::new(),
            driver_gene_data: Vec::new(),
            gene_copy_numbers: Vec::new(),
            driver_output: Vec::new(),
            sample_ploidy: 0.0,
            is_male: false,
        }
    }

    /// setSamplePurityData
    pub fn set_sample_purity_data(&mut self, ploidy: f64, is_male: bool) {
        self.sample_ploidy = ploidy;
        self.is_male = is_male;
    }
    /// addDriverGene (tests): a driver and its gene copy number in place of the PURPLE files
    pub fn add_driver_gene(&mut self, driver: DriverCatalog, gene_cn: GeneCopyNumber) {
        self.driver_catalog.clear();
        self.driver_gene_data.clear();
        self.gene_copy_numbers.clear();
        self.driver_catalog.push(driver);
        self.gene_copy_numbers.push(gene_cn);
    }
    /// clearResults
    pub fn clear_results(&mut self, clear_drivers: bool) {
        self.driver_gene_data.clear();
        if clear_drivers {
            self.driver_catalog.clear();
        }
    }

    fn load_data_from_file(&mut self, purple_dir: &str, sample_id: &str) {
        self.driver_catalog.clear();
        self.driver_gene_data.clear();
        self.gene_copy_numbers.clear();
        let r = (|| -> anyhow::Result<()> {
            let pc = read_purity_context(purple_dir, sample_id)?;
            self.sample_ploidy = pc.ploidy;
            self.is_male = pc.gender.starts_with("MALE");
            let cat = read_driver_catalog(&driver_catalog_somatic_filename(purple_dir.trim_end_matches('/'), sample_id))?;
            const PURPLE_SOMATIC: [&str; 6] = ["AMP", "PARTIAL_AMP", "DEL", "HET_DEL", "LOH", "MUTATION"];
            self.driver_catalog.extend(cat.into_iter().filter(|x| PURPLE_SOMATIC.contains(&x.driver.as_str()) && x.reported_status == "REPORTED"));
            self.gene_copy_numbers.extend(read_gene_copy_numbers(&gene_copy_number_filename(purple_dir.trim_end_matches('/'), sample_id))?);
            Ok(())
        })();
        if let Err(e) = r {
            hmf_common::error!("failed to load driver catalog or purity context: {}", e);
        }
    }

    fn find_gene_copy_number(&self, gene_name: &str) -> Option<&GeneCopyNumber> {
        self.gene_copy_numbers.iter().find(|g| g.gene_name == gene_name)
    }

    /// annotateSVs
    pub fn annotate_svs(&mut self, sample_id: &str, chr_breakend_map: &JavaHashMap<String, Vec<Be>>, cn: &CnDataLoader, vis: &mut VisSampleData) {
        self.driver_output.clear();
        if !self.purple_data_path.is_empty() {
            let dir = self.purple_data_path.replace('*', sample_id);
            self.load_data_from_file(&dir, sample_id);
        }
        for i in 0..self.driver_catalog.len() {
            let dc = self.driver_catalog[i].clone();
            if !driver_type_handled(&dc) {
                continue;
            }
            let Some(gene) = self.gene_cache.gene_data_by_name(&dc.gene) else {
                hmf_common::warn!("driver gene({}) data not found", dc.gene);
                continue;
            };
            let Some(trans) = self.gene_cache.transcript_data(&gene.gene_id, &dc.transcript) else {
                hmf_common::warn!("driver gene({}) trans data not found", dc.gene);
                continue;
            };
            let region = match self.find_gene_copy_number(&dc.gene) {
                Some(g) => GeneCopyNumberRegion::from_gene_cn(g),
                None => match calc_gene_copy_number_region(trans, cn.chr_cn_data_map.get(&gene.chromosome)) {
                    Some(r) => r,
                    None => {
                        hmf_common::warn!("sample({}) gene({}) min copy number data not found for driver", sample_id, dc.gene);
                        continue;
                    }
                },
            };
            let mut dg = DriverGeneData::new(dc.clone(), gene, trans, region);
            if dc.driver == "DEL" {
                self.annotate_delete_event(&mut dg, cn);
            } else if dc.category == "TSG" && dc.biallelic {
                annotate_biallelic_event(&mut dg, cn);
            } else if dc.driver == "AMP" || dc.driver == "PARTIAL_AMP" {
                let list = chr_breakend_map.get(&dg.gene.chromosome);
                self.annotate_amplification(&mut dg, list, cn);
            }
            self.write_driver_data(&dg, vis);
            self.driver_gene_data.push(dg);
        }
        let dis_del = self.find_disruptive_del_drivers(chr_breakend_map, cn);
        for dg in &dis_del {
            self.write_driver_data(dg, vis);
        }
        self.driver_gene_data.extend(dis_del);
    }

    fn write_driver_data(&mut self, dg: &DriverGeneData, vis: &mut VisSampleData) {
        if !dg.trans.is_canonical {
            return;
        }
        for e in &dg.events {
            let cluster_id = e.cluster().map(|c| c.id()).unwrap_or(-1);
            self.driver_output.push(LinxDriver { cluster_id, gene: dg.gene.gene_name.clone(), event_type: e.event_type });
            vis.add_gene_exon_data(VisGeneData::new(
                cluster_id,
                &dg.gene.gene_id,
                &dg.gene.gene_name,
                "",
                0,
                &dg.gene.chromosome,
                VisGeneAnnotationType::DRIVER,
            ));
        }
    }

    // DeletionDrivers

    fn annotate_delete_event(&self, dg: &mut DriverGeneData, cn: &CnDataLoader) {
        if dg.cn_region.is_germline || dg.cn_region.min_copy_number >= 0.51 {
            return;
        }
        let min_start = dg.cn_region.region_start;
        let min_end = dg.cn_region.region_end;
        if CnDataLoader::expect_single_chromosome(self.is_male, &dg.gene.chromosome) {
            for hl in &cn.hom_loss_data {
                if hl.pos_start <= min_end && hl.pos_end >= min_start {
                    let mut e = DriverGeneEvent::new(ET::DEL);
                    e.add_sv_breakend_pair(hl.breakend(true), hl.breakend(false), "DEL");
                    dg.events.push(e);
                    break;
                }
            }
            return;
        }
        for loh in &cn.loh_event_data {
            if loh.chromosome != dg.gene.chromosome {
                continue;
            }
            if loh.pos_start > min_end || loh.pos_end < min_start {
                continue;
            }
            if (loh.pos_start < min_end && loh.pos_end > min_start) || !loh.double_sv_event() {
                add_loh_event(dg, loh);
            } else {
                let mut e = DriverGeneEvent::new(ET::DEL);
                e.add_sv_breakend_pair(loh.breakend(true), loh.breakend(false), "DEL");
                dg.events.push(e);
            }
            for hl in loh.hom_loss_events() {
                if hl.pos_start <= min_end && hl.pos_end >= min_start {
                    let mut e = DriverGeneEvent::new(ET::DEL);
                    e.add_sv_breakend_pair(hl.breakend(true), hl.breakend(false), "DEL");
                    dg.events.push(e);
                    break;
                }
            }
            break;
        }
    }

    fn find_disruptive_del_drivers(&mut self, chr_breakend_map: &JavaHashMap<String, Vec<Be>>, cn: &CnDataLoader) -> Vec<DriverGeneData> {
        let mut del_gene_ids: Vec<String> =
            self.driver_gene_data.iter().filter(|x| x.driver.driver == "DEL").map(|x| x.gene.gene_id.clone()).collect();
        let mut out = Vec::new();
        for (_, list) in chr_breakend_map.iter() {
            for breakend in list {
                let del_type = breakend.forward_orient() && breakend.db_link().is_some();
                let dup_type = breakend.reverse_orient() && breakend.sv_type() == DUP;
                if !del_type && !dup_type {
                    continue;
                }
                let genes: Vec<BreakendGene> = breakend
                    .sv()
                    .genes_list(breakend.uses_start())
                    .into_iter()
                    .filter(|x| !del_gene_ids.iter().any(|g| g == x.gene_id()))
                    .filter(|x| self.reportable_disruption_gene_transcripts.contains_key(x.gene_id()))
                    .collect();
                for gene in &genes {
                    for t in gene.transcripts() {
                        if !t.is_disruptive() {
                            continue;
                        }
                        if !t.is_canonical() && !self.reportable_disruption_gene_transcripts[gene.gene_id()].iter().any(|x| x == t.trans_name()) {
                            continue;
                        }
                        if del_type {
                            self.check_del_disruption(breakend, &t, &mut del_gene_ids, &mut out, cn);
                        } else {
                            self.check_dup_disruption(breakend, &t, &mut del_gene_ids, &mut out, cn);
                        }
                    }
                }
            }
        }
        out
    }

    fn check_del_disruption(&mut self, breakend: &Be, t: &BreakendTrans, del_gene_ids: &mut Vec<String>, out: &mut Vec<DriverGeneData>, cn: &CnDataLoader) {
        let db = breakend.db_link().expect("deletion bridge");
        let mut cn_low = breakend.copy_number_low_side();
        let other = db.other_breakend(breakend);
        if db.length() < 0 {
            cn_low -= other.jcn();
        }
        if cn_low >= TOTAL_CN_LOSS {
            return;
        }
        let gene_id = t.gene().gene_id();
        let Some(matching) = other.sv().genes_list(other.uses_start()).into_iter().find(|x| x.gene_id() == gene_id) else { return };
        if !matching.transcripts().iter().any(|x| x.trans_name() == t.trans_name() && x.is_disruptive()) {
            return;
        }
        del_gene_ids.push(gene_id.to_string());
        if let Some(mut dg) = self.create_hom_disruption_driver(t.gene(), t.trans_data, "HOM_DEL_DISRUPTION", cn) {
            let mut e = DriverGeneEvent::new(ET::HOM_DEL_DISRUPTION);
            e.add_sv_breakend_pair(Some(breakend.clone()), Some(other.clone()), "DB");
            e.cluster = breakend.cluster();
            dg.events.push(e);
            out.push(dg);
        }
    }

    fn check_dup_disruption(&mut self, breakend: &Be, t: &BreakendTrans, del_gene_ids: &mut Vec<String>, out: &mut Vec<DriverGeneData>, cn: &CnDataLoader) {
        let gene_id = t.gene().gene_id();
        if !breakend.sv().genes_list(!breakend.uses_start()).iter().any(|x| x.gene_id() == gene_id) {
            return;
        }
        let other = breakend.other_breakend().expect("DUP other breakend");
        let Some(matching) = other.sv().genes_list(other.uses_start()).into_iter().find(|x| x.gene_id() == gene_id) else { return };
        if !matching.transcripts().iter().any(|x| x.trans_name() == t.trans_name() && x.is_disruptive()) {
            return;
        }
        if is_homozygous_dup_disruption(breakend, &other, breakend.jcn()) {
            del_gene_ids.push(gene_id.to_string());
            if let Some(mut dg) = self.create_hom_disruption_driver(t.gene(), t.trans_data, "HOM_DUP_DISRUPTION", cn) {
                let mut e = DriverGeneEvent::new(ET::HOM_DUP_DISRUPTION);
                e.add_sv_breakend_pair(Some(breakend.clone()), Some(other.clone()), "DUP");
                e.cluster = breakend.cluster();
                dg.events.push(e);
                out.push(dg);
            }
        }
    }

    fn create_hom_disruption_driver(&mut self, gene: &BreakendGene, trans: Trans, driver_type: &str, cn: &CnDataLoader) -> Option<DriverGeneData> {
        let gcn = self.find_gene_copy_number(gene.gene_name()).cloned();
        let dg = self.driver_genes.get(&gene.gene_name().to_string())?;
        let record = DriverCatalog {
            chromosome: gene.chromosome(),
            chromosome_band: gene.gene_data.karyotype_band.clone(),
            gene: gene.gene_name().to_string(),
            transcript: trans.trans_name.clone(),
            is_canonical: trans.is_canonical,
            driver: driver_type.to_string(),
            category: dg.likelihood_type.name().to_string(),
            likelihood_method: "DISRUPTION".into(),
            reported_status: if dg.report_disruption { "REPORTED".into() } else { "NOT_REPORTED".into() },
            driver_likelihood: 1.0,
            missense: 0,
            nonsense: 0,
            splice: 0,
            inframe: 0,
            frameshift: 0,
            biallelic: true,
            min_copy_number: gcn.as_ref().map(|g| g.min_copy_number).unwrap_or(0.0),
            max_copy_number: gcn.as_ref().map(|g| g.max_copy_number).unwrap_or(0.0),
        };
        self.driver_catalog.push(record.clone());
        let Some(gene_data) = self.gene_cache.gene_data_by_name(gene.gene_name()) else {
            hmf_common::warn!("gene({}) no Ensembl gene data found", gene.gene_name());
            return None;
        };
        let region = match &gcn {
            Some(g) => GeneCopyNumberRegion::from_gene_cn(g),
            // Java goes on with a null region where none is found
            None => calc_gene_copy_number_region(trans, cn.chr_cn_data_map.get(&gene_data.chromosome)).expect("gene copy number region"),
        };
        Some(DriverGeneData::new(record, gene_data, trans, region))
    }

    // AmplificationDrivers

    fn annotate_amplification(&self, dg: &mut DriverGeneData, breakends: Option<&Vec<Be>>, cn: &CnDataLoader) {
        self.check_chromosome_amplification(dg, cn);
        check_cluster_amplification(dg, breakends);
        if dg.events.is_empty() {
            dg.events.push(DriverGeneEvent::new(ET::GAIN));
        }
    }

    fn check_chromosome_amplification(&self, dg: &mut DriverGeneData, cn: &CnDataLoader) {
        let Some(tc) = cn.chr_ends_cn_map.get(&dg.gene.chromosome) else { return };
        let telomere_cn = tc.telomere_p_arm.min(tc.telomere_q_arm);
        let centromere_cn = tc.centromere_p_arm.min(tc.centromere_q_arm);
        let chromosome_cn = centromere_cn.min(telomere_cn);
        if !has_short_arm(&dg.gene.chromosome) && chromosome_cn / self.sample_ploidy > 2.0 {
            let mut e = DriverGeneEvent::new(ET::GAIN_CHR);
            e.copy_number_gain = chromosome_cn - self.sample_ploidy;
            dg.events.push(e);
        }
        let centromere_change = if dg.arm == Arm::P {
            tc.centromere_p_arm.min(tc.telomere_p_arm) - tc.centromere_q_arm
        } else {
            tc.centromere_q_arm.min(tc.telomere_q_arm) - tc.centromere_p_arm
        };
        if centromere_change > 0.0 && !copy_numbers_equal(centromere_change, 0.0) {
            let mut e = DriverGeneEvent::new(ET::GAIN_ARM);
            e.copy_number_gain = centromere_change;
            dg.events.push(e);
        }
    }
}

/// HumanChromosome.hasShortArm: the acrocentric chromosomes. Java compares the short names ("13") with the
/// chromosome as given, which GRCh38's "chr13" never equals (LX5).
pub(crate) fn has_short_arm(chromosome: &str) -> bool {
    let chr = if crate::compat::java("LX5") { chromosome } else { hmf_common::regions::strip_chr_prefix(chromosome) };
    matches!(chr, "13" | "14" | "15" | "21" | "22")
}

fn driver_type_handled(d: &DriverCatalog) -> bool {
    if d.driver == "DEL" || d.driver == "AMP" || d.driver == "PARTIAL_AMP" {
        return true;
    }
    d.driver == "MUTATION" && d.category == "TSG" && d.biallelic
}

fn add_loh_event(dg: &mut DriverGeneData, loh: &crate::cn::Loh) {
    if loh.double_sv_event() {
        let mut e = DriverGeneEvent::new(ET::LOH);
        e.add_sv_breakend_pair(loh.breakend(true), loh.breakend(false), "DEL");
        dg.events.push(e);
    } else if loh.arm_loss() {
        dg.events.push(DriverGeneEvent::new(ET::LOH_ARM));
    } else if loh.chromosome_loss() {
        dg.events.push(DriverGeneEvent::new(ET::LOH_CHR));
    } else if loh.is_sv_event_any() {
        let t = if loh.telomere_loss() { ET::LOH_SV_TELO } else { ET::LOH_SV_CENTRO };
        let mut e = DriverGeneEvent::new(t);
        e.add_sv_breakend_pair(loh.breakend(true), loh.breakend(false), "ARM_SV");
        dg.events.push(e);
    }
}

fn annotate_biallelic_event(dg: &mut DriverGeneData, cn: &CnDataLoader) {
    let (Some(cs), Some(ce)) = (dg.trans.coding_start, dg.trans.coding_end) else { return };
    for loh in &cn.loh_event_data {
        if loh.chromosome != dg.gene.chromosome {
            continue;
        }
        if loh.pos_start > ce || loh.pos_end < cs {
            continue;
        }
        add_loh_event(dg, loh);
        break;
    }
}

pub fn is_homozygous_dup_disruption(breakend: &Be, other: &Be, jcn: f64) -> bool {
    let s = breakend.copy_number_low_side();
    let e = other.copy_number_low_side();
    let threshold = (jcn * (1.0 + MAX_COPY_NUM_DIFF_PERC)).max(jcn + MAX_COPY_NUM_DIFF);
    s < threshold && e < threshold
}

struct OpposingSegment {
    breakends: Vec<Be>,
    remaining_cn_change: f64,
}

struct DriverAmpData {
    traverse_up: bool,
    breakend_count: i32,
    segment_count: i32,
    start_copy_number: f64,
    net_cn_change: f64,
}

/// returns the index of the opposing segment, or None where it adds no opposing copy number
fn find_or_create_opposing_segment(segs: &mut Vec<OpposingSegment>, list: &[Be], start: &Be, traverse_up: bool, stop: i32) -> Option<usize> {
    if let Some(i) = segs.iter().position(|x| contains(&x.breakends, start)) {
        return Some(i);
    }
    let start_cn = start.get_copy_number(traverse_up);
    let mut end_cn = start.get_copy_number(!traverse_up);
    let mut index = start.chr_pos_index();
    let mut seg_breakends = vec![start.clone()];
    loop {
        index += if traverse_up { 1 } else { -1 };
        if index < 0 || index as usize >= list.len() {
            break;
        }
        let b = &list[index as usize];
        if (traverse_up && b.position() > stop) || (!traverse_up && b.position() < stop) {
            break;
        }
        if b.cluster() == start.cluster() {
            seg_breakends.push(b.clone());
            end_cn = b.get_copy_number(!traverse_up);
        } else {
            break;
        }
    }
    let net = end_cn - start_cn;
    if net > 0.0 || copy_numbers_equal(net, 0.0) {
        segs.push(OpposingSegment { breakends: seg_breakends, remaining_cn_change: 0.0 });
        return None;
    }
    segs.push(OpposingSegment { breakends: seg_breakends, remaining_cn_change: -net });
    Some(segs.len() - 1)
}

fn check_cluster_for_amplification(
    dg: &DriverGeneData, list: &[Be], start: &Be, traverse_up: bool, segs: &mut Vec<OpposingSegment>,
) -> Option<DriverAmpData> {
    let trans_start = dg.trans.trans_start;
    let trans_end = dg.trans.trans_end;
    let start_cn = start.get_copy_number(traverse_up);
    let target = start.cluster();
    let mut in_segment = true;
    let mut seg_start_cn = start_cn;
    let mut net = 0.0;
    let segment_count = 0;
    let mut breakend_count = 0;
    let mut index = start.chr_pos_index();
    let mut breakend: Option<Be> = None;
    let mut seg_start = start.clone();
    loop {
        if breakend.is_some() {
            index += if traverse_up { 1 } else { -1 };
        }
        if index < 0 || index as usize >= list.len() {
            break;
        }
        let b = list[index as usize].clone();
        breakend = Some(b.clone());
        if (traverse_up && b.position() > trans_start) || (!traverse_up && b.position() < trans_end) {
            break;
        }
        let cluster = b.cluster();
        if cluster != target {
            if in_segment {
                let end_cn = if b.arm() == seg_start.arm() {
                    b.get_copy_number(traverse_up)
                } else {
                    let prev = (index + if traverse_up { -1 } else { 1 }) as usize;
                    list[prev].get_copy_number(!traverse_up)
                };
                net += end_cn - seg_start_cn;
                in_segment = false;
            }
            let opp = find_or_create_opposing_segment(segs, list, &b, traverse_up, if traverse_up { trans_start } else { trans_end });
            if let Some(oi) = opp {
                let o = &mut segs[oi];
                if net > 0.0 && o.remaining_cn_change > 0.0 {
                    if o.remaining_cn_change > net {
                        o.remaining_cn_change = (o.remaining_cn_change - net).max(0.0);
                        net = 0.0;
                    } else {
                        net -= o.remaining_cn_change;
                        o.remaining_cn_change = 0.0;
                    }
                }
            }
        } else {
            breakend_count += 1;
            if !in_segment {
                seg_start_cn = b.get_copy_number(traverse_up);
                seg_start = b.clone();
                in_segment = true;
            }
        }
    }
    let mut cluster_change = net;
    if in_segment {
        cluster_change += dg.cn_region.min_copy_number - seg_start_cn;
    }
    if cluster_change > 0.0 && !copy_numbers_equal(cluster_change, 0.0) && !copy_numbers_equal(start_cn + cluster_change, start_cn) {
        return Some(DriverAmpData {
            traverse_up,
            breakend_count,
            segment_count,
            start_copy_number: start_cn,
            net_cn_change: cluster_change,
        });
    }
    None
}

fn check_cluster_amplification(dg: &mut DriverGeneData, list: Option<&Vec<Be>>) {
    let Some(list) = list else { return };
    if list.is_empty() {
        return;
    }
    let trans_start = dg.trans.trans_start;
    let trans_end = dg.trans.trans_end;
    let start_index: i32 = 0;
    let end_index: i32 = list.len() as i32 - 1;
    // a HashMap keyed by cluster (identity) in Java (LX1)
    let mut amp: IdMap<Cluster, DriverAmpData> = IdMap::new();
    for i in 0..=1 {
        let traverse_up = i == 0;
        let mut processed: Vec<Cluster> = Vec::new();
        let mut segs: Vec<OpposingSegment> = Vec::new();
        let mut index = if traverse_up { start_index } else { end_index };
        let mut breakend: Option<Be> = None;
        loop {
            if breakend.is_some() {
                index += if traverse_up { 1 } else { -1 };
            }
            if index < start_index || index > end_index {
                break;
            }
            let b = list[index as usize].clone();
            breakend = Some(b.clone());
            if (traverse_up && b.position() > trans_end) || (!traverse_up && b.position() < trans_start) {
                break;
            }
            let cluster = b.cl();
            if cluster.has_linking_line_elements() || contains(&processed, &cluster) {
                continue;
            }
            if cluster.sv_count() == 1 && cluster.sv0().sv_type() == DEL {
                continue;
            }
            processed.push(cluster.clone());
            let Some(a) = check_cluster_for_amplification(dg, list, &b, traverse_up, &mut segs) else { continue };
            let replace = amp.get(&cluster).is_none_or(|e| e.net_cn_change < a.net_cn_change);
            if replace {
                amp.insert(cluster, a);
            }
        }
    }
    if amp.is_empty() {
        return;
    }
    let mut max_cn_change = amp.iter().map(|(_, a)| a.net_cn_change).fold(f64::NEG_INFINITY, f64::max);
    let mut index = 0;
    while index < dg.events.len() {
        let gain = dg.events[index].copy_number_gain;
        if gain < MIN_AMP_PERCENT_VS_MAX * max_cn_change {
            dg.events.remove(index);
        } else {
            max_cn_change = max_cn_change.max(gain);
            index += 1;
        }
    }
    for (cluster, a) in amp.iter() {
        if a.net_cn_change / max_cn_change < MIN_AMP_PERCENT_VS_MAX {
            continue;
        }
        let mut e = DriverGeneEvent::new(ET::GAIN);
        e.cluster = Some(cluster.clone());
        e.copy_number_gain = a.net_cn_change;
        e.sv_info = format!(
            "{};{};{};{}",
            a.traverse_up,
            a.breakend_count,
            a.segment_count,
            hmf_common::fmt::java_fixed(a.start_copy_number, 1)
        );
        dg.events.push(e);
    }
}
