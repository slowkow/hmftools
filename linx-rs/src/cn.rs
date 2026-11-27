//! Package `cn`: PURPLE's copy number segments as LINX holds them (SvCNData), the loader that links them to
//! the SVs and finds LOH and homozygous-loss events (CnDataLoader), the junction copy number recalculation
//! (CnJcnCalcs) and the copy number segments made from the SVs alone (CnSegmentBuilder, germline mode).

use std::cell::{Cell, RefCell};
use std::rc::Rc;
use std::sync::Mutex;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::{self, PurityContext};
use hmf_common::sv::StructuralVariantData;
use hmf_common::sv::StructuralVariantType::*;

use crate::types::*;
use crate::utils::{jmax, jmin, jround};

/// SvCNData: a copy number segment
pub struct SvCNData {
    id: i32,
    pub chromosome: String,
    pub start_pos: i32,
    pub end_pos: i32,
    pub copy_number: f64,
    pub seg_start: String,
    pub seg_end: String,
    pub method: String,
    pub baf_count: i32,
    pub observed_baf: f64,
    pub actual_baf: f64,
    pub depth_window_count: i32,
    index: Cell<i32>,
    sv_data: RefCell<Option<Rc<StructuralVariantData>>>,
    sv_link_on_start: Cell<bool>,
}

handle!(SvCnData, SvCNData);

impl SvCNData {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: i32, chromosome: &str, start_pos: i32, end_pos: i32, copy_number: f64, seg_start: &str, seg_end: &str, baf_count: i32,
        actual_baf: f64, depth_window_count: i32,
    ) -> SvCNData {
        SvCNData {
            id,
            chromosome: chromosome.to_string(),
            start_pos,
            end_pos,
            copy_number,
            seg_start: seg_start.to_string(),
            seg_end: seg_end.to_string(),
            method: "STRUCTURAL_VARIANT".to_string(),
            baf_count,
            observed_baf: 0.0,
            actual_baf,
            depth_window_count,
            index: Cell::new(0),
            sv_data: RefCell::new(None),
            sv_link_on_start: Cell::new(false),
        }
    }

    pub fn from_record(r: &purple_files::PurpleCopyNumber, id: i32) -> SvCNData {
        SvCNData {
            id,
            chromosome: r.chromosome.clone(),
            start_pos: r.start,
            end_pos: r.end,
            copy_number: r.average_tumor_copy_number,
            seg_start: r.segment_start_support.clone(),
            seg_end: r.segment_end_support.clone(),
            method: r.method.clone(),
            baf_count: r.baf_count,
            observed_baf: r.average_observed_baf,
            actual_baf: r.average_actual_baf,
            depth_window_count: r.depth_window_count,
            index: Cell::new(0),
            sv_data: RefCell::new(None),
            sv_link_on_start: Cell::new(false),
        }
    }

    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn position(&self, use_start: bool) -> i32 {
        if use_start { self.start_pos } else { self.end_pos }
    }
    pub fn length(&self) -> i32 {
        self.end_pos - self.start_pos
    }
    pub fn index(&self) -> i32 {
        self.index.get()
    }
    pub fn set_index(&self, i: i32) {
        self.index.set(i)
    }
    pub fn major_allele_jcn(&self) -> f64 {
        self.actual_baf * self.copy_number
    }
    pub fn minor_allele_jcn(&self) -> f64 {
        jmax((1.0 - self.actual_baf) * self.copy_number, 0.0)
    }
    pub fn set_structural_variant_data(&self, sv: &Rc<StructuralVariantData>, link_on_start: bool) {
        *self.sv_data.borrow_mut() = Some(sv.clone());
        self.sv_link_on_start.set(link_on_start);
    }
    pub fn structural_variant_data(&self) -> Option<Rc<StructuralVariantData>> {
        self.sv_data.borrow().clone()
    }
    pub fn sv_link_on_start(&self) -> bool {
        self.sv_link_on_start.get()
    }
    pub fn matches_segment(&self, segment: &str, is_start: bool) -> bool {
        if is_start { self.seg_start == segment } else { self.seg_end == segment }
    }
    /// matchesSV: the segment boundary is an SV type (or MULTIPLE)
    pub fn matches_sv(&self, is_start: bool) -> bool {
        is_sv_segment(if is_start { &self.seg_start } else { &self.seg_end })
    }
}

/// SvCNData.isSvSegment
pub fn is_sv_segment(segment: &str) -> bool {
    matches!(segment, "BND" | "INV" | "INS" | "DEL" | "DUP" | "SGL" | "MULTIPLE")
}

impl std::fmt::Display for SvCNData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "chr({}) pos({} - {}) segs({} - {}) cn({}) baf({})",
            self.chromosome,
            self.start_pos,
            self.end_pos,
            self.seg_start,
            self.seg_end,
            crate::utils::ff(self.copy_number, 2),
            crate::utils::ff(self.actual_baf, 2)
        )
    }
}

pub const TELOMERE: &str = "TELOMERE";
pub const CENTROMERE: &str = "CENTROMERE";
pub const SEG_NONE: &str = "NONE";
pub const SEG_UNKNOWN: &str = "UNKNOWN";

pub const CN_DATA_NO_SV: i32 = -1;

/// LohEvent
pub struct LohEvent {
    pub chromosome: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub seg_start: String,
    pub seg_end: String,
    pub seg_count: i32,
    pub start_sv: i32,
    pub end_sv: i32,
    hom_loss_events: RefCell<Vec<HomLoss>>,
    start_cn_data: RefCell<Option<SvCnData>>,
    end_cn_data: RefCell<Option<SvCnData>>,
    is_valid: Cell<bool>,
    breakend_start: RefCell<Option<Be>>,
    breakend_end: RefCell<Option<Be>>,
}

handle!(Loh, LohEvent);

impl LohEvent {
    #[allow(clippy::too_many_arguments)]
    pub fn new(chr: &str, pos_start: i32, pos_end: i32, seg_start: &str, seg_end: &str, seg_count: i32, start_sv: i32, end_sv: i32) -> LohEvent {
        LohEvent {
            chromosome: chr.to_string(),
            pos_start,
            pos_end,
            seg_start: seg_start.to_string(),
            seg_end: seg_end.to_string(),
            seg_count,
            start_sv,
            end_sv,
            hom_loss_events: RefCell::new(Vec::new()),
            start_cn_data: RefCell::new(None),
            end_cn_data: RefCell::new(None),
            is_valid: Cell::new(true),
            breakend_start: RefCell::new(None),
            breakend_end: RefCell::new(None),
        }
    }
    pub fn set_breakend(&self, breakend: Option<Be>, is_start: bool) {
        if is_start {
            *self.breakend_start.borrow_mut() = breakend;
        } else {
            *self.breakend_end.borrow_mut() = breakend;
        }
    }
    pub fn breakend(&self, is_start: bool) -> Option<Be> {
        if is_start { self.breakend_start.borrow().clone() } else { self.breakend_end.borrow().clone() }
    }
    pub fn set_cn_data(&self, start: Option<SvCnData>, end: Option<SvCnData>) {
        *self.start_cn_data.borrow_mut() = start;
        *self.end_cn_data.borrow_mut() = end;
    }
    pub fn cn_data(&self, is_start: bool) -> Option<SvCnData> {
        if is_start { self.start_cn_data.borrow().clone() } else { self.end_cn_data.borrow().clone() }
    }
    pub fn matched_both_svs(&self) -> bool {
        self.breakend_start.borrow().is_some() && self.breakend_end.borrow().is_some()
    }
    pub fn same_sv(&self) -> bool {
        match (&*self.breakend_start.borrow(), &*self.breakend_end.borrow()) {
            (Some(s), Some(e)) => s.sv() == e.sv(),
            // Java: mBreakendEnd.getSV() throws if the end is unset
            (Some(_), None) => panic!("LohEvent.sameSV without an end breakend"),
            _ => false,
        }
    }
    pub fn length(&self) -> i32 {
        self.pos_end - self.pos_start
    }
    pub fn add_hom_loss_events(&self, events: &[HomLoss]) {
        self.hom_loss_events.borrow_mut().extend(events.iter().cloned());
        let v = !self.has_incomplete_hom_loss_events();
        self.is_valid.set(v);
    }
    pub fn hom_loss_events(&self) -> Vec<HomLoss> {
        self.hom_loss_events.borrow().clone()
    }
    pub fn has_incomplete_hom_loss_events(&self) -> bool {
        self.hom_loss_events.borrow().iter().any(|x| !x.matched_both_svs() || !x.same_sv())
    }
    pub fn clustered(&self) -> bool {
        self.matched_both_svs() && self.breakend(true).unwrap().cluster() == self.breakend(false).unwrap().cluster()
    }
    pub fn is_valid(&self) -> bool {
        self.is_valid.get()
    }
    pub fn set_is_valid(&self, t: bool) {
        self.is_valid.set(t)
    }
    pub fn arm_loss(&self) -> bool {
        self.arm_loss_of(Arm::P) || self.arm_loss_of(Arm::Q)
    }
    pub fn arm_loss_of(&self, arm: Arm) -> bool {
        (arm == Arm::P && self.seg_start == TELOMERE && self.seg_end == CENTROMERE)
            || (arm == Arm::Q && self.seg_start == CENTROMERE && self.seg_end == TELOMERE)
    }
    pub fn telomere_loss(&self) -> bool {
        self.seg_start == TELOMERE || self.seg_end == TELOMERE
    }
    pub fn centromere_loss(&self) -> bool {
        self.seg_start == CENTROMERE || self.seg_end == CENTROMERE
    }
    pub fn chromosome_loss(&self) -> bool {
        self.seg_start == TELOMERE && self.seg_end == TELOMERE
    }
    pub fn whole_arm_loss(&self) -> bool {
        self.arm_loss() || self.chromosome_loss()
    }
    pub fn is_sv_event(&self, use_start: bool) -> bool {
        if use_start { self.start_sv != CN_DATA_NO_SV } else { self.end_sv != CN_DATA_NO_SV }
    }
    pub fn double_sv_event(&self) -> bool {
        self.is_sv_event(true) && self.is_sv_event(false)
    }
    /// isSvEvent(): either end at an SV
    pub fn is_sv_event_any(&self) -> bool {
        self.is_sv_event(true) || self.is_sv_event(false)
    }
}

impl std::fmt::Display for LohEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.breakend(true).map(|b| b.sv().id().to_string()).unwrap_or_else(|| "none".into());
        let e = self.breakend(false).map(|b| b.sv().id().to_string()).unwrap_or_else(|| "none".into());
        write!(f, "chr({}) segs({} -> {}) pos({} -> {}) SVs({} & {})", self.chromosome, self.seg_start, self.seg_end, self.pos_start, self.pos_end, s, e)
    }
}
impl std::fmt::Display for Loh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// HomLossEvent
pub struct HomLossEvent {
    pub chromosome: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub seg_start: String,
    pub seg_end: String,
    pub start_sv: i32,
    pub end_sv: i32,
    breakend_start: RefCell<Option<Be>>,
    breakend_end: RefCell<Option<Be>>,
}

handle!(HomLoss, HomLossEvent);

impl HomLossEvent {
    pub fn new(chr: &str, pos_start: i32, pos_end: i32, seg_start: &str, seg_end: &str, start_sv: i32, end_sv: i32) -> HomLossEvent {
        HomLossEvent {
            chromosome: chr.to_string(),
            pos_start,
            pos_end,
            seg_start: seg_start.to_string(),
            seg_end: seg_end.to_string(),
            start_sv,
            end_sv,
            breakend_start: RefCell::new(None),
            breakend_end: RefCell::new(None),
        }
    }
    pub fn set_breakend(&self, breakend: Option<Be>, is_start: bool) {
        if is_start {
            *self.breakend_start.borrow_mut() = breakend;
        } else {
            *self.breakend_end.borrow_mut() = breakend;
        }
    }
    pub fn clear_breakends(&self) {
        *self.breakend_start.borrow_mut() = None;
        *self.breakend_end.borrow_mut() = None;
    }
    pub fn breakend(&self, is_start: bool) -> Option<Be> {
        if is_start { self.breakend_start.borrow().clone() } else { self.breakend_end.borrow().clone() }
    }
    pub fn matched_both_svs(&self) -> bool {
        self.breakend_start.borrow().is_some() && self.breakend_end.borrow().is_some()
    }
    pub fn same_sv(&self) -> bool {
        match (&*self.breakend_start.borrow(), &*self.breakend_end.borrow()) {
            (Some(s), Some(e)) => s.sv() == e.sv(),
            (Some(_), None) => panic!("HomLossEvent.sameSV without an end breakend"),
            _ => false,
        }
    }
    pub fn clustered(&self) -> bool {
        self.matched_both_svs() && self.breakend(true).unwrap().cluster() == self.breakend(false).unwrap().cluster()
    }
    pub fn cluster(&self) -> Option<Cluster> {
        if self.matched_both_svs() { self.breakend(true).unwrap().cluster() } else { None }
    }
}

impl std::fmt::Display for HomLossEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.breakend(true).map(|b| b.sv().id().to_string()).unwrap_or_else(|| "none".into());
        let e = self.breakend(false).map(|b| b.sv().id().to_string()).unwrap_or_else(|| "none".into());
        write!(f, "chr({}) segs({} -> {}) pos({} -> {}) SVs({} & {})", self.chromosome, self.seg_start, self.seg_end, self.pos_start, self.pos_end, s, e)
    }
}

/// JcnCalcData
#[derive(Clone, Copy, Debug)]
pub struct JcnCalcData {
    pub jcn_estimate: f64,
    pub jcn_uncertainty: f64,
    pub valid: bool,
}

/// TelomereCentromereCnData
#[derive(Clone, Copy, Debug)]
pub struct TelomereCentromereCnData {
    pub telomere_p_arm: f64,
    pub telomere_q_arm: f64,
    pub centromere_p_arm: f64,
    pub centromere_q_arm: f64,
}

pub const MIN_LOH_CN: f64 = 0.5;
pub const TOTAL_CN_LOSS: f64 = 0.5;
const REMOTE_SV_DISTANCE: i32 = 1000000;

/// CnDataLoader (with CnJcnCalcs, whose maps it shares)
pub struct CnDataLoader {
    purple_data_path: String,
    record_id: i32,
    sv_data_list: Vec<Rc<StructuralVariantData>>,
    pub chr_cn_data_map: JavaHashMap<String, Vec<SvCnData>>,
    pub sv_id_cn_data_map: JavaHashMap<i32, [Option<SvCnData>; 2]>,
    pub sv_jcn_calc_map: JavaHashMap<i32, JcnCalcData>,
    pub purity_context: Option<PurityContext>,
    pub loh_event_data: Vec<Loh>,
    pub hom_loss_data: Vec<HomLoss>,
    pub chr_ends_cn_map: JavaHashMap<String, TelomereCentromereCnData>,
}

impl CnDataLoader {
    pub fn new(purple_data_path: &str) -> CnDataLoader {
        let p = if !purple_data_path.is_empty() {
            if purple_data_path.ends_with('/') { purple_data_path.to_string() } else { format!("{purple_data_path}/") }
        } else {
            String::new()
        };
        CnDataLoader {
            purple_data_path: p,
            record_id: 0,
            sv_data_list: Vec::new(),
            chr_cn_data_map: JavaHashMap::new(),
            sv_id_cn_data_map: JavaHashMap::new(),
            sv_jcn_calc_map: JavaHashMap::new(),
            purity_context: None,
            loh_event_data: Vec::new(),
            hom_loss_data: Vec::new(),
            chr_ends_cn_map: JavaHashMap::new(),
        }
    }

    pub fn load_sample_data(&mut self, sample_id: &str, sv_records: &[Rc<StructuralVariantData>]) {
        self.sv_data_list = sv_records.to_vec();
        self.load_copy_number_data(sample_id);
        self.create_chr_copy_number_map();
        self.link_copy_number_and_sv_data(sample_id);
        self.find_loh_events(sample_id);
        self.calculate_adjusted_jcn(sample_id);
    }

    fn load_copy_number_data(&mut self, sample_id: &str) {
        self.chr_cn_data_map.clear();
        let mut records = Vec::new();
        if !self.purple_data_path.is_empty() {
            let dir = if self.purple_data_path.contains('*') {
                self.purple_data_path.replace('*', sample_id)
            } else {
                self.purple_data_path.clone()
            };
            let r = purple_files::read_purple_copy_numbers(&purple_files::purple_copy_number_filename(&dir, sample_id))
                .and_then(|cn| purple_files::read_purity_context(&dir, sample_id).map(|p| (cn, p)));
            match r {
                Ok((cn, p)) => {
                    records = cn;
                    self.purity_context = Some(p);
                }
                Err(e) => {
                    hmf_common::error!("failed to load purity context: {:#}", e);
                    return;
                }
            }
        }
        let mut current_chromosome = String::new();
        for rec in &records {
            self.record_id += 1;
            let cn_data = SvCnData::new(SvCNData::from_record(rec, self.record_id));
            if current_chromosome != cn_data.chromosome {
                current_chromosome = cn_data.chromosome.clone();
                self.chr_cn_data_map.insert(current_chromosome.clone(), Vec::new());
            }
            let list = self.chr_cn_data_map.get_mut(&current_chromosome).unwrap();
            cn_data.set_index(list.len() as i32);
            list.push(cn_data);
        }
    }

    fn link_copy_number_and_sv_data(&mut self, sample_id: &str) {
        self.sv_id_cn_data_map.clear();
        let mut unmatched_svs = 0;
        for sv_data in &self.sv_data_list {
            if sv_data.filter == hmf_common::sv::PON_FILTER_PON {
                continue;
            }
            for be in [SE_START, SE_END] {
                let is_start = is_start(be);
                let sv_chromosome = if is_start { &sv_data.start_chromosome } else { &sv_data.end_chromosome };
                let sv_position = if is_start { sv_data.start_position } else { sv_data.end_position };
                let sv_orientation = if is_start { sv_data.start_orientation } else { sv_data.end_orientation };
                let Some(cn_list) = self.chr_cn_data_map.get(sv_chromosome) else { continue };
                if cn_list.is_empty() {
                    continue;
                }
                let mut found = false;
                for cn_data in cn_list {
                    let mut cn_position = cn_data.start_pos;
                    if sv_orientation == 1 {
                        cn_position -= 1;
                    }
                    if sv_position != cn_position {
                        continue;
                    }
                    if cn_data.structural_variant_data().is_none() {
                        cn_data.set_structural_variant_data(sv_data, is_start);
                    }
                    if !self.sv_id_cn_data_map.contains_key(&sv_data.id) {
                        self.sv_id_cn_data_map.insert(sv_data.id, [None, None]);
                    }
                    self.sv_id_cn_data_map.get_mut(&sv_data.id).unwrap()[be] = Some(cn_data.clone());
                    found = true;
                }
                if !found && sv_data.sv_type() != INS {
                    hmf_common::debug!(
                        "SV({} chr={} pos={} type={}) {} unmatched)",
                        sv_data.id,
                        sv_chromosome,
                        sv_position,
                        sv_data.sv_type(),
                        if is_start { "start" } else { "end" }
                    );
                    unmatched_svs += 1;
                }
            }
        }
        if unmatched_svs > 0 {
            hmf_common::warn!("sample({}) has {} unmatched CN-SV segments", sample_id, unmatched_svs);
        }
    }

    fn is_single_variant(cn_data: &SvCNData) -> bool {
        let Some(sv) = cn_data.structural_variant_data() else { return false };
        let sv_pos_start = if sv.start_orientation == -1 { cn_data.start_pos } else { cn_data.start_pos - 1 };
        let sv_pos_end = cn_data.end_pos + 1;
        sv.start_position == sv_pos_start && sv.end_position == sv_pos_end
    }

    pub fn is_male_sample(&self) -> bool {
        self.purity_context.as_ref().is_some_and(|p| p.gender.starts_with("MALE"))
    }

    /// CnDataLoader.expectSingleChromosome: X and Y in a male sample. Java compares with "X" and "Y", which
    /// GRCh38's "chrX" and "chrY" never equal (LX2).
    pub fn expect_single_chromosome(is_male: bool, chromosome: &str) -> bool {
        let chr = if crate::compat::java("LX2") { chromosome } else { hmf_common::regions::strip_chr_prefix(chromosome) };
        is_male && (chr == "X" || chr == "Y")
    }

    fn find_loh_events(&mut self, sample_id: &str) {
        self.loh_event_data.clear();
        self.hom_loss_data.clear();
        let mut is_loh_section = false;
        let mut loh_min_cn = 0.0;
        let mut last_min_cn = 0.0;
        let mut _prior_cn = 0.0;
        let mut loh_segments = 0;
        let mut loh_start_cn_data: Option<SvCnData> = None;
        let mut loh_section_count = 0;
        let mut loh_svs_matched_count = 0;
        let mut loh_on_start_telomere = false;
        let mut total_loss = false;
        let mut loh_hom_loss_events: Vec<HomLoss> = Vec::new();

        let entries: Vec<(String, Vec<SvCnData>)> = self.chr_cn_data_map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (chromosome, cn_list) in &entries {
            let mut index = 0usize;
            while index < cn_list.len() {
                let cn_data = &cn_list[index];
                let min_cn = cn_data.minor_allele_jcn();
                let new_chromosome = index == 0;
                let mut reset = new_chromosome;

                if cn_data.copy_number < TOTAL_CN_LOSS {
                    let next = if index < cn_list.len() - 1 { Some(&cn_list[index + 1]) } else { None };
                    let sv = find_sv_data(cn_data, 1);
                    let next_sv = next.and_then(|n| find_sv_data(n, -1));
                    let hom_loss = HomLoss::new(HomLossEvent::new(
                        chromosome,
                        cn_data.start_pos,
                        cn_data.end_pos,
                        &cn_data.seg_start,
                        &cn_data.seg_end,
                        sv.map(|s| s.id).unwrap_or(CN_DATA_NO_SV),
                        next_sv.map(|s| s.id).unwrap_or(CN_DATA_NO_SV),
                    ));
                    self.hom_loss_data.push(hom_loss.clone());
                    if is_loh_section || loh_on_start_telomere {
                        loh_hom_loss_events.push(hom_loss);
                    }
                }

                if is_loh_section || loh_on_start_telomere {
                    let mut loh_regained = min_cn >= MIN_LOH_CN;
                    if loh_regained && cn_data.matches_sv(true) && find_sv_data(cn_data, -1).is_none() {
                        loh_regained = false;
                        hmf_common::debug!("skipping LOH end segment({}) with no SV match", cn_data.0);
                    }
                    if loh_regained || reset {
                        if loh_regained && cn_data.end_pos - cn_data.start_pos <= SHORT_TI_LENGTH && index < cn_list.len() - 1 {
                            let next = &cn_list[index + 1];
                            let next_min_cn = next.minor_allele_jcn();
                            if next.end_pos - cn_data.start_pos > REMOTE_SV_DISTANCE
                                && next_min_cn < MIN_LOH_CN
                                && loh_start_cn_data.as_ref().is_some_and(|s| cn_data.start_pos - s.start_pos > REMOTE_SV_DISTANCE)
                            {
                                // a 'continue' in Java: lastMinCN is not updated
                                index += 1;
                                continue;
                            }
                        }
                        if loh_on_start_telomere || total_loss {
                            self.process_loh_data(chromosome, loh_start_cn_data.clone(), Some(cn_data.clone()), loh_segments, false, &loh_hom_loss_events);
                            loh_on_start_telomere = false;
                            total_loss = false;
                        } else {
                            loh_svs_matched_count += self.process_loh_data(
                                chromosome,
                                loh_start_cn_data.clone(),
                                Some(cn_data.clone()),
                                loh_segments,
                                false,
                                &loh_hom_loss_events,
                            );
                            loh_section_count += 1;
                        }
                        reset = true;
                    } else if cn_data.matches_segment(TELOMERE, false) {
                        self.process_loh_data(chromosome, loh_start_cn_data.clone(), Some(cn_data.clone()), loh_segments, true, &loh_hom_loss_events);
                        reset = true;
                    } else if cn_data.copy_number < TOTAL_CN_LOSS {
                        if cn_data.matches_segment("DEL", true) && cn_data.matches_segment("DEL", false) && Self::is_single_variant(cn_data) {
                            // a total CN loss matching a single SV
                        } else {
                            total_loss = true;
                        }
                    } else {
                        loh_segments += 1;
                        loh_min_cn = jmin(loh_min_cn, min_cn);
                    }
                } else if !reset {
                    let mut loh_lost = min_cn < MIN_LOH_CN;
                    let sv = find_sv_data(cn_data, 1);
                    if loh_lost && cn_data.matches_sv(true) && sv.is_none() {
                        loh_lost = false;
                    }
                    if loh_lost {
                        is_loh_section = true;
                        loh_segments = 1;
                        loh_start_cn_data = Some(cn_data.clone());
                        loh_min_cn = min_cn;
                        _prior_cn = last_min_cn;
                        if !loh_on_start_telomere && cn_data.matches_segment(TELOMERE, false) {
                            self.process_loh_data(chromosome, loh_start_cn_data.clone(), loh_start_cn_data.clone(), loh_segments, true, &[]);
                            reset = true;
                        }
                    }
                }

                if reset {
                    is_loh_section = false;
                    loh_on_start_telomere = false;
                    loh_segments = 0;
                    loh_start_cn_data = None;
                    total_loss = false;
                    loh_hom_loss_events.clear();
                    if new_chromosome && min_cn < MIN_LOH_CN {
                        loh_on_start_telomere = true;
                        loh_start_cn_data = Some(cn_data.clone());
                    }
                }
                last_min_cn = min_cn;
                index += 1;
            }
        }
        let _ = loh_min_cn;
        hmf_common::debug!("sample({}) LOH sections({}) fullMatched({})", sample_id, loh_section_count, loh_svs_matched_count);
    }

    #[allow(clippy::too_many_arguments)]
    fn process_loh_data(
        &mut self, chromosome: &str, start_data: Option<SvCnData>, end_data: Option<SvCnData>, seg_count: i32, incomplete: bool,
        loh_hom_loss_events: &[HomLoss],
    ) -> i32 {
        if Self::expect_single_chromosome(self.is_male_sample(), chromosome) {
            return 0;
        }
        // Java dereferences the start segment; a null one throws
        let start = start_data.expect("LOH start segment");
        let mut end = end_data;
        let start_sv = find_sv_data(&start, 1);
        let end_sv: Option<Rc<StructuralVariantData>> = if incomplete {
            None
        } else if !start.matches_segment(TELOMERE, false) {
            find_sv_data(end.as_ref().expect("LOH end segment"), -1)
        } else {
            end = Some(start.clone());
            None
        };
        let end = end.expect("LOH end segment");

        if let (Some(ss), Some(es)) = (&start_sv, &end_sv) {
            if ss.id != es.id {
                for (data, sv) in [(&start, ss), (&end, es)] {
                    if data.matches_segment("DEL", true) || data.matches_segment("DUP", true) {
                        if let Some(items) = self.sv_id_cn_data_map.get(&sv.id) {
                            if let (Some(s), Some(e)) = (&items[SE_START], &items[SE_END]) {
                                if s.index() == e.index() - 1 {
                                    hmf_common::debug!(
                                        "segs start({}) and end({}) skipped since bounded by simpleSV({})",
                                        start.0,
                                        end.0,
                                        sv.id
                                    );
                                    return 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        let loh = Loh::new(LohEvent::new(
            chromosome,
            start.start_pos,
            if incomplete { end.end_pos } else { end.start_pos },
            &start.seg_start,
            if incomplete { &end.seg_end } else { &end.seg_start },
            seg_count,
            start_sv.as_ref().map(|s| s.id).unwrap_or(CN_DATA_NO_SV),
            end_sv.as_ref().map(|s| s.id).unwrap_or(CN_DATA_NO_SV),
        ));
        loh.add_hom_loss_events(loh_hom_loss_events);
        loh.set_cn_data(Some(start), Some(end));
        self.loh_event_data.push(loh);
        if start_sv.is_some() && end_sv.is_some() { 1 } else { 0 }
    }

    pub fn create_chr_copy_number_map(&mut self) {
        self.chr_ends_cn_map.clear();
        let mut m = JavaHashMap::new();
        for (chr, list) in self.chr_cn_data_map.iter() {
            let mut d = TelomereCentromereCnData { telomere_p_arm: 0.0, telomere_q_arm: 0.0, centromere_p_arm: 0.0, centromere_q_arm: 0.0 };
            for cn in list {
                if cn.seg_start == TELOMERE {
                    d.telomere_p_arm = cn.copy_number;
                }
                if cn.seg_end == CENTROMERE {
                    d.centromere_p_arm = cn.copy_number;
                }
                if cn.seg_start == CENTROMERE {
                    d.centromere_q_arm = cn.copy_number;
                }
                if cn.seg_end == TELOMERE {
                    d.telomere_q_arm = cn.copy_number;
                }
            }
            m.insert(chr.clone(), d);
        }
        self.chr_ends_cn_map = m;
    }

    fn cn_segment(&self, chromosome: &str, index: i32) -> Option<SvCnData> {
        let list = self.chr_cn_data_map.get(&chromosome.to_string())?;
        if index < 0 || index as usize >= list.len() {
            return None;
        }
        Some(list[index as usize].clone())
    }

    /// CnJcnCalcs.calculateAdjustedJcn
    pub fn calculate_adjusted_jcn(&mut self, sample_id: &str) {
        self.sv_jcn_calc_map.clear();
        let mut results: Vec<(i32, JcnCalcData)> = Vec::new();
        for (sv_id, pair) in self.sv_id_cn_data_map.iter() {
            let sv_id = *sv_id;
            let Some(cn_start) = pair[SE_START].clone().filter(|c| c.structural_variant_data().is_some()) else {
                hmf_common::error!("SV({}) missing start copy number data or unlinked SV data", sv_id);
                continue;
            };
            let mut sv = cn_start.structural_variant_data().unwrap();
            if sv.id != sv_id {
                if let Some(s) = self.sv_data_list.iter().find(|s| s.id == sv_id) {
                    sv = s.clone();
                }
            }
            // Java dereferences the previous segment unchecked
            let cn_start_prev = self.cn_segment(&cn_start.chromosome, cn_start.index() - 1).expect("previous CN segment");
            let cn_end_next = pair[SE_END].clone();
            let cn_end = cn_end_next.as_ref().and_then(|n| self.cn_segment(&n.chromosome, n.index() - 1));
            let tumor_read_count_start = sv.start_tumor_variant_fragment_count;
            let max_cn_start = jmax(cn_start_prev.copy_number, cn_start.copy_number);
            let max_cn_end = match (&cn_end, &cn_end_next) {
                (Some(e), Some(n)) => jmax(e.copy_number, n.copy_number),
                _ => 0.0,
            };
            let mut start_depth = [cn_start_prev.depth_window_count, cn_start.depth_window_count];
            let mut cn_chg_start = sv.adjusted_start_copy_number_change;
            let mut end_depth =
                [cn_end.as_ref().map(|c| c.depth_window_count).unwrap_or(0), cn_end_next.as_ref().map(|c| c.depth_window_count).unwrap_or(0)];
            let mut cn_chg_end = sv.adjusted_end_copy_number_change;
            if let Some(e) = &cn_end {
                if sv.sv_type() == INV && cn_start.index() == e.index() {
                    let n = cn_end_next.as_ref().unwrap();
                    if cn_start_prev.depth_window_count > cn_start.depth_window_count && n.depth_window_count > cn_start.depth_window_count {
                        let cn_change = (cn_start_prev.copy_number - n.copy_number).abs() * 0.5;
                        end_depth[0] = end_depth[1];
                        start_depth[1] = start_depth[0];
                        cn_chg_start = cn_change;
                        cn_chg_end = cn_change;
                    }
                }
            }
            let calc = calc_adjusted_jcn_values(
                cn_chg_start,
                cn_chg_end,
                tumor_read_count_start,
                sv.junction_copy_number,
                max_cn_start,
                max_cn_end,
                &start_depth,
                if cn_end.is_some() { Some(&end_depth) } else { None },
            );
            if !calc.valid {
                hmf_common::debug!(
                    "sample({}) svID({} type={}) unexpected JCN(est={} unc={})",
                    sample_id,
                    sv.id,
                    sv.sv_type(),
                    crate::utils::dstr(calc.jcn_estimate),
                    crate::utils::dstr(calc.jcn_uncertainty)
                );
            }
            results.push((sv.id, calc));
        }
        for (id, c) in results {
            self.sv_jcn_calc_map.insert(id, c);
        }
    }
}

/// CnDataLoader.findSvData: the SV starting at this segment with the required orientation
fn find_sv_data(cn_data: &SvCNData, required_orient: i8) -> Option<Rc<StructuralVariantData>> {
    if cn_data.matches_segment(SEG_UNKNOWN, true) || cn_data.matches_segment(TELOMERE, true) || cn_data.matches_segment(CENTROMERE, true) {
        return None;
    }
    let sv_position = if required_orient == -1 { cn_data.start_pos } else { cn_data.start_pos - 1 };
    let sv = cn_data.structural_variant_data()?;
    if cn_data.sv_link_on_start() && sv.start_orientation == required_orient && sv.start_position == sv_position {
        Some(sv)
    } else if !cn_data.sv_link_on_start() && sv.end_orientation == required_orient && sv.end_position == sv_position {
        Some(sv)
    } else {
        None
    }
}

const ABS_UNCERTAINTY: f64 = 0.15;
const RELATIVE_UNCERTAINTY: f64 = 0.10;
const ADDITIONAL_ABS_UNCERTAINTY: f64 = 0.4;
const ADDITIONAL_REL_UNCERTAINTY: f64 = 0.15;
const PROPORTION_CNCHANGE_USED_IN_JCN_UNC: f64 = 0.5;
pub const READ_COUNT_PROB_LOW: f64 = 0.005;
pub const READ_COUNT_PROB_HIGH: f64 = 0.995;
const READ_COUNT_PROB_MAX: i32 = 2000;

/// CnJcnCalcs.calcAdjustedJcnValues
#[allow(clippy::too_many_arguments)]
pub fn calc_adjusted_jcn_values(
    cn_chg_start: f64, cn_chg_end: f64, tumor_read_count: i32, jcn: f64, max_cn_start: f64, max_cn_end: f64, start_depth: &[i32; 2],
    end_depth: Option<&[i32; 2]>,
) -> JcnCalcData {
    let cn_unc_start = calc_copy_number_side_uncertainty(max_cn_start, start_depth);
    let cn_unc_end = end_depth.map(|d| calc_copy_number_side_uncertainty(max_cn_end, d)).unwrap_or(0.0);
    let mut observations: Vec<f64> = Vec::new();
    let mut uncertainties: Vec<f64> = Vec::new();
    if cn_unc_start > 0.0 {
        observations.push(cn_chg_start);
        uncertainties.push(cn_unc_start);
    }
    if cn_unc_end > 0.0 {
        observations.push(cn_chg_end);
        uncertainties.push(cn_unc_end);
    }
    if tumor_read_count > 0 {
        let low = calc_poisson_observed_given_prob_cached(tumor_read_count, true) as f64;
        let high = calc_poisson_observed_given_prob_cached(tumor_read_count, false) as f64;
        let rc_adjusted_jcn = jcn * ((low + high) * 0.5) / tumor_read_count as f64;
        let mut jcn_uncertainty = rc_adjusted_jcn * (high - low) / tumor_read_count as f64 * 0.5;
        let cn_unc_factor = if !uncertainties.is_empty() {
            hmf_common::java_math::double_stream_sum(uncertainties.iter().copied()) / uncertainties.len() as f64
        } else {
            0.0
        };
        jcn_uncertainty += cn_unc_factor * PROPORTION_CNCHANGE_USED_IN_JCN_UNC;
        if jcn_uncertainty > 0.0 {
            observations.push(rc_adjusted_jcn);
            uncertainties.push(jcn_uncertainty);
        }
    }
    let (est_jcn, est_uncertainty);
    if observations.len() == 1 {
        est_jcn = observations[0];
        est_uncertainty = uncertainties[0];
    } else {
        let mut sum_uncertainty = 0.0;
        let mut sum_observed_uncertainty = 0.0;
        for i in 0..observations.len() {
            let inv = 1.0 / uncertainties[i].powi(2);
            sum_uncertainty += inv;
            sum_observed_uncertainty += observations[i] * inv;
        }
        est_jcn = sum_observed_uncertainty / sum_uncertainty;
        let mut adj_uncertainty = 0.0;
        for i in 0..observations.len() {
            let inv = 1.0 / uncertainties[i].powi(2);
            adj_uncertainty += inv * jmax(observations[i] - est_jcn, uncertainties[i] / 2.0).powi(2);
        }
        est_uncertainty = (observations.len() as f64 / (observations.len() as f64 - 1.0) * adj_uncertainty / sum_uncertainty).sqrt();
    }
    if est_jcn.is_nan() || est_jcn <= 0.0 || est_uncertainty.is_nan() || est_uncertainty <= 0.0 {
        return JcnCalcData { jcn_estimate: 0.0, jcn_uncertainty: 0.0, valid: false };
    }
    JcnCalcData { jcn_estimate: est_jcn, jcn_uncertainty: est_uncertainty, valid: true }
}

fn calc_copy_number_side_uncertainty(copy_number: f64, depth: &[i32; 2]) -> f64 {
    let min_depth_count = jmax(depth[0].min(depth[1]) as f64, 0.1);
    let mut u = jmax(copy_number * RELATIVE_UNCERTAINTY, ABS_UNCERTAINTY);
    u += jmax(ADDITIONAL_ABS_UNCERTAINTY, ADDITIONAL_REL_UNCERTAINTY * copy_number) / min_depth_count.sqrt();
    u
}

/// The read count cache of CnJcnCalcs (establishReadCountCache), filled as values are asked for: the same
/// numbers as Java's eager cache, since both come from calcPoissonObservedGivenProb
static READ_COUNT_CACHE: Mutex<Vec<Option<[i32; 2]>>> = Mutex::new(Vec::new());

fn calc_poisson_observed_given_prob_cached(expected: i32, use_low: bool) -> i32 {
    if expected <= 0 {
        return 0;
    }
    if expected <= READ_COUNT_PROB_MAX {
        let mut cache = READ_COUNT_CACHE.lock().unwrap();
        if cache.is_empty() {
            cache.resize(READ_COUNT_PROB_MAX as usize + 1, None);
        }
        let e = cache[expected as usize].get_or_insert_with(|| {
            [calc_poisson_observed_given_prob(expected, READ_COUNT_PROB_LOW), calc_poisson_observed_given_prob(expected, READ_COUNT_PROB_HIGH)]
        });
        return if use_low { e[0] } else { e[1] };
    }
    calc_poisson_observed_given_prob(expected, if use_low { READ_COUNT_PROB_LOW } else { READ_COUNT_PROB_HIGH })
}

/// CnJcnCalcs.calcPoissonObservedGivenProb: the count with this cumulative probability, by bisection
pub fn calc_poisson_observed_given_prob(expected: i32, required_prob: f64) -> i32 {
    if expected <= 0 {
        return 0;
    }
    let poisson = hmf_common::poisson::Poisson::new(expected as f64);
    let max_iterations = 20;
    let mut iterations = 0;
    let ref_count = 25.0;
    let ref_range_perc = 0.44;
    let range_perc = ref_range_perc / (expected as f64 / ref_count).sqrt();
    let range = range_perc * expected as f64;
    let (mut upper, mut lower): (i32, i32);
    if required_prob > 0.5 {
        upper = jround(expected as f64 + range * 1.5) as i32;
        lower = jround(jmax(expected as f64 + range * 0.2, 0.0)) as i32;
    } else {
        upper = jround(expected as f64 - range * 0.2) as i32;
        lower = jround(jmax(expected as f64 - range * 1.2, 0.0)) as i32;
    }
    let mut test_value = ((lower + upper) as f64 * 0.5) as i32;
    let mut current_prob = poisson.cumulative_probability(test_value);
    let mut prob_diff = 0.0;
    while iterations < max_iterations {
        prob_diff = (required_prob - current_prob).abs() / required_prob;
        if prob_diff < 0.001 {
            break;
        }
        if current_prob > required_prob {
            if test_value <= lower + 1 {
                break;
            }
            upper = test_value;
            test_value = jround((test_value + lower) as f64 * 0.5) as i32;
        } else {
            if test_value >= upper - 1 {
                break;
            }
            lower = test_value;
            test_value = jround((test_value + upper) as f64 * 0.5) as i32;
        }
        current_prob = poisson.cumulative_probability(test_value);
        iterations += 1;
    }
    if iterations >= max_iterations {
        hmf_common::warn!(
            "max iterations reached: value({}) test({}) prob({} diff={})",
            expected,
            test_value,
            crate::utils::ff(current_prob, 4),
            crate::utils::ff(prob_diff, 4)
        );
    }
    test_value
}

/// SvCNData.setSvCopyNumberData: the JCN recalculation and the segments either side of each breakend
pub fn set_sv_copy_number_data(
    sv_list: &[Sv], jcn_calc: Option<&JavaHashMap<i32, JcnCalcData>>, sv_id_cn: &JavaHashMap<i32, [Option<SvCnData>; 2]>,
    chr_cn: &JavaHashMap<String, Vec<SvCnData>>,
) {
    if jcn_calc.is_none_or(|m| m.is_empty()) && sv_id_cn.is_empty() {
        return;
    }
    let mut current_chromosome = String::new();
    let mut cn_list: Option<&Vec<SvCnData>> = None;
    for var in sv_list {
        if let Some(m) = jcn_calc {
            if let Some(j) = m.get(&var.id()) {
                var.set_jcn_recalc_data(j.jcn_estimate - j.jcn_uncertainty, j.jcn_estimate + j.jcn_uncertainty);
            }
        }
        let Some(pair) = sv_id_cn.get(&var.id()) else { continue };
        for be in [SE_START, SE_END] {
            if var.is_sgl_breakend() && be == SE_END {
                continue;
            }
            let is_start = is_start(be);
            if current_chromosome != var.chromosome(is_start) {
                current_chromosome = var.chromosome(is_start).to_string();
                cn_list = chr_cn.get(&current_chromosome);
            }
            let Some(post) = &pair[be] else { continue };
            let Some(list) = cn_list else { continue };
            // Java's List.get throws for index -1
            let prev = list[(post.index() - 1) as usize].clone();
            var.set_copy_number_segments(is_start, Some(prev), Some(post.clone()));
        }
    }
}

/// CnSegmentBuilder
pub struct CnSegmentBuilder {
    other_allele_jcn: f64,
    undisrupted_allele_jcn: f64,
}

impl Default for CnSegmentBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl CnSegmentBuilder {
    pub fn new() -> Self {
        CnSegmentBuilder { other_allele_jcn: 1.0, undisrupted_allele_jcn: 0.0 }
    }
    pub fn set_allele_ploidies(&mut self, other: f64, undisrupted: f64) {
        self.other_allele_jcn = other;
        self.undisrupted_allele_jcn = undisrupted;
    }

    fn calc_actual_baf(&self, copy_number: f64) -> f64 {
        if copy_number == 0.0 {
            return 0.0;
        }
        let b = jmax(copy_number - self.other_allele_jcn, 0.0);
        if b >= self.other_allele_jcn { b / copy_number } else { self.other_allele_jcn / copy_number }
    }

    /// createCopyNumberData: segments from the SVs' JCNs (used by the tests)
    pub fn create_copy_number_data(&self, loader: &mut CnDataLoader, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
        loader.chr_cn_data_map.clear();
        loader.sv_id_cn_data_map.clear();
        let mut cn_id = 0;
        for (chromosome, breakend_list) in chr_breakend_map.iter() {
            let mut cn_list: Vec<SvCnData> = Vec::new();
            let net_sv_jcn = jmax(
                hmf_common::java_math::double_stream_sum(breakend_list.iter().map(|x| x.jcn() * x.orientation() as f64)),
                0.0,
            );
            let mut current = self.other_allele_jcn + self.undisrupted_allele_jcn + net_sv_jcn;
            let centromere_position = crate::utils::get_chromosomal_arm_length(chromosome, Arm::P);
            let chromosome_length = crate::utils::get_chromosome_length(chromosome);
            let push = |list: &mut Vec<SvCnData>, d: SvCNData| -> SvCnData {
                let d = SvCnData::new(d);
                d.set_index(list.len() as i32);
                list.push(d.clone());
                d
            };
            for i in 0..breakend_list.len() {
                let breakend = &breakend_list[i];
                let var = breakend.sv();
                let sv_data = var.sv_data_rc().clone();
                let jcn = var.jcn();
                let jcn_change = -jcn * breakend.orientation() as f64;
                let tstr = var.sv_type().name();
                if i == 0 {
                    if breakend.sv_type() == DUP && breakend_list.get(i + 1).is_some_and(|n| n.sv() == breakend.sv()) {
                        current += jcn_change;
                    } else {
                        current += jmax(-jcn_change, 0.0);
                    }
                    if current < 0.0 {
                        hmf_common::error!("invalid copy number({}) at telomere", crate::utils::dstr(current));
                        return;
                    }
                    let baf = self.calc_actual_baf(current);
                    if breakend.arm() == Arm::Q {
                        push(&mut cn_list, SvCNData::new(cn_id, chromosome, 0, centromere_position, current, TELOMERE, CENTROMERE, 1, baf, 100));
                        cn_id += 1;
                        push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, centromere_position, breakend.position() - 1, current, CENTROMERE, tstr, 1, baf, 100),
                        );
                        cn_id += 1;
                    } else {
                        push(&mut cn_list, SvCNData::new(cn_id, chromosome, 0, breakend.position() - 1, current, TELOMERE, tstr, 1, baf, 100));
                        cn_id += 1;
                    }
                }
                current += jcn_change;
                if current < 0.0 {
                    hmf_common::warn!(
                        "invalid copy number({}) from jcnChange({}) at breakend({})",
                        crate::utils::dstr(current),
                        crate::utils::format_jcn(jcn_change),
                        breakend
                    );
                    current = 0.0;
                }
                let baf = self.calc_actual_baf(current);
                let cn_data;
                if i < breakend_list.len() - 1 {
                    let next = &breakend_list[i + 1];
                    if breakend.arm() == Arm::P && next.arm() == Arm::Q {
                        cn_data = push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, breakend.position(), centromere_position - 1, current, tstr, CENTROMERE, 1, baf, 100),
                        );
                        cn_id += 1;
                        cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                        push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, centromere_position, next.position() - 1, current, CENTROMERE, next.sv_type().name(), 1, baf, 100),
                        );
                        cn_id += 1;
                    } else {
                        cn_data = push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, breakend.position(), next.position() - 1, current, tstr, next.sv_type().name(), 1, baf, 100),
                        );
                        cn_id += 1;
                        cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                    }
                } else if breakend.arm() == Arm::P {
                    cn_data = push(
                        &mut cn_list,
                        SvCNData::new(cn_id, chromosome, breakend.position(), centromere_position - 1, current, tstr, CENTROMERE, 1, baf, 100),
                    );
                    cn_id += 1;
                    cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                    push(&mut cn_list, SvCNData::new(cn_id, chromosome, centromere_position, chromosome_length, current, CENTROMERE, TELOMERE, 1, 0.5, 100));
                    cn_id += 1;
                } else {
                    cn_data = push(
                        &mut cn_list,
                        SvCNData::new(cn_id, chromosome, breakend.position(), chromosome_length, current, tstr, TELOMERE, 1, baf, 100),
                    );
                    cn_id += 1;
                    cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                }
                if !loader.sv_id_cn_data_map.contains_key(&var.id()) {
                    loader.sv_id_cn_data_map.insert(var.id(), [None, None]);
                }
                loader.sv_id_cn_data_map.get_mut(&var.id()).unwrap()[if breakend.uses_start() { SE_START } else { SE_END }] = Some(cn_data);
                let be_cn = if breakend.forward_orient() { current + jcn } else { current };
                var.set_copy_number_data(breakend.uses_start(), be_cn, jcn);
            }
            loader.chr_cn_data_map.insert(chromosome.clone(), cn_list);
        }
    }

    /// createGermlineCopyNumberData: a segment for each breakend, irrespective of the breakends around it
    pub fn create_germline_copy_number_data(&self, loader: &mut CnDataLoader, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
        loader.chr_cn_data_map.clear();
        loader.sv_id_cn_data_map.clear();
        let default_cn = 2.0;
        let default_baf = 0.5;
        let mut cn_id = 0;
        for (chromosome, breakend_list) in chr_breakend_map.iter() {
            let mut cn_list: Vec<SvCnData> = Vec::new();
            let centromere_position = crate::utils::get_chromosomal_arm_length(chromosome, Arm::P);
            let chromosome_length = crate::utils::get_chromosome_length(chromosome);
            let push = |list: &mut Vec<SvCnData>, d: SvCNData| -> SvCnData {
                let d = SvCnData::new(d);
                d.set_index(list.len() as i32);
                list.push(d.clone());
                d
            };
            for i in 0..breakend_list.len() {
                let breakend = &breakend_list[i];
                let var = breakend.sv();
                let sv_data = var.sv_data_rc().clone();
                let variant_cn = breakend.copy_number();
                let variant_baf = if variant_cn >= 2.0 { (variant_cn - 1.0) / variant_cn } else { jmin(1.0, 1.0 / variant_cn) };
                let tstr = var.sv_type().name();
                if i == 0 {
                    if breakend.arm() == Arm::Q {
                        push(&mut cn_list, SvCNData::new(cn_id, chromosome, 0, centromere_position, default_cn, TELOMERE, CENTROMERE, 1, default_baf, 100));
                        cn_id += 1;
                        push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, centromere_position, breakend.position() - 1, default_cn, CENTROMERE, tstr, 1, default_baf, 100),
                        );
                        cn_id += 1;
                    } else {
                        push(&mut cn_list, SvCNData::new(cn_id, chromosome, 0, breakend.position() - 1, default_cn, TELOMERE, tstr, 1, default_baf, 100));
                        cn_id += 1;
                    }
                }
                let cn_data;
                if i < breakend_list.len() - 1 {
                    let next = &breakend_list[i + 1];
                    if breakend.arm() == Arm::P && next.arm() == Arm::Q {
                        cn_data = push(
                            &mut cn_list,
                            SvCNData::new(cn_id, chromosome, breakend.position(), centromere_position - 1, variant_cn, tstr, CENTROMERE, 1, variant_baf, 100),
                        );
                        cn_id += 1;
                        cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                        push(
                            &mut cn_list,
                            SvCNData::new(
                                cn_id,
                                chromosome,
                                centromere_position,
                                next.position() - 1,
                                default_cn,
                                CENTROMERE,
                                next.sv_type().name(),
                                1,
                                default_baf,
                                100,
                            ),
                        );
                        cn_id += 1;
                    } else {
                        cn_data = push(
                            &mut cn_list,
                            SvCNData::new(
                                cn_id,
                                chromosome,
                                breakend.position(),
                                next.position() - 1,
                                variant_cn,
                                tstr,
                                next.sv_type().name(),
                                1,
                                variant_baf,
                                100,
                            ),
                        );
                        cn_id += 1;
                        cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                    }
                } else if breakend.arm() == Arm::P {
                    cn_data = push(
                        &mut cn_list,
                        SvCNData::new(cn_id, chromosome, breakend.position(), centromere_position - 1, variant_cn, tstr, CENTROMERE, 1, variant_baf, 100),
                    );
                    cn_id += 1;
                    cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                    push(
                        &mut cn_list,
                        SvCNData::new(cn_id, chromosome, centromere_position, chromosome_length, default_cn, CENTROMERE, TELOMERE, 1, default_baf, 100),
                    );
                    cn_id += 1;
                } else {
                    cn_data = push(
                        &mut cn_list,
                        SvCNData::new(cn_id, chromosome, breakend.position(), chromosome_length, variant_cn, tstr, TELOMERE, 1, variant_baf, 100),
                    );
                    cn_id += 1;
                    cn_data.set_structural_variant_data(&sv_data, breakend.uses_start());
                }
                if !loader.sv_id_cn_data_map.contains_key(&var.id()) {
                    loader.sv_id_cn_data_map.insert(var.id(), [None, None]);
                }
                loader.sv_id_cn_data_map.get_mut(&var.id()).unwrap()[if breakend.uses_start() { SE_START } else { SE_END }] = Some(cn_data);
                let d = var.sv_data();
                if breakend.uses_start() {
                    var.set_copy_number_data(true, d.adjusted_start_copy_number, d.adjusted_start_copy_number_change);
                } else {
                    var.set_copy_number_data(false, d.adjusted_end_copy_number, d.adjusted_end_copy_number_change);
                }
            }
            loader.chr_cn_data_map.insert(chromosome.clone(), cn_list);
        }
    }

    /// setSamplePurity: a purity context of the given purity, ploidy and gender
    pub fn set_sample_purity(loader: &mut CnDataLoader, purity: f64, ploidy: f64, gender: &str) {
        loader.purity_context = Some(PurityContext { purity, ploidy, gender: gender.to_string(), whole_genome_duplication: false });
    }
}
