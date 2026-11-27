//! The LINX data model (package `types`, and SvChain from `chaining`): SVs, their breakends, linked pairs
//! (templated insertions), deletion bridges, clusters, chains, arm groups and arm clusters.
//!
//! Java's object graph is kept as it is: each object lives behind an `Rc` (the handle types `Sv`, `Be`,
//! `Pair`, `Db`, `Cluster`, `Chain`), compared by identity as Java's `==` and `List.contains` compare
//! them, with its mutable fields in `Cell`s and `RefCell`s. Lists handed out by the getters are copies
//! (of the handles) unless the method says otherwise, so that no borrow is held across a call.

use std::cell::{Cell, Ref, RefCell, RefMut};
use std::rc::Rc;

use hmf_common::java_hash::JavaHashMap;

use hmf_common::sv::{StructuralVariantData, StructuralVariantType};
use hmf_common::sv::StructuralVariantType::*;

use crate::cn::SvCnData;
use crate::gene::BreakendGene;

pub const SE_START: usize = 0;
pub const SE_END: usize = 1;

#[inline]
pub fn se(is_start: bool) -> usize {
    if is_start { SE_START } else { SE_END }
}
#[inline]
pub fn is_start(se: usize) -> bool {
    se == SE_START
}

/// A handle type: an `Rc` compared, hashed and keyed by identity
macro_rules! handle {
    ($name:ident, $inner:ty) => {
        #[derive(Clone)]
        pub struct $name(pub std::rc::Rc<$inner>);
        impl PartialEq for $name {
            fn eq(&self, o: &Self) -> bool {
                std::rc::Rc::ptr_eq(&self.0, &o.0)
            }
        }
        impl Eq for $name {}
        impl std::hash::Hash for $name {
            fn hash<H: std::hash::Hasher>(&self, h: &mut H) {
                (std::rc::Rc::as_ptr(&self.0) as *const u8 as usize).hash(h)
            }
        }
        impl std::ops::Deref for $name {
            type Target = $inner;
            fn deref(&self) -> &$inner {
                &self.0
            }
        }
        impl hmf_common::java_id_hash::IdKey for $name {
            fn id_addr(&self) -> usize {
                std::rc::Rc::as_ptr(&self.0) as *const u8 as usize
            }
        }
        impl $name {
            pub fn new(v: $inner) -> Self {
                $name(std::rc::Rc::new(v))
            }
        }
    };
}
pub(crate) use handle;

handle!(Sv, SvVarData);
handle!(Be, SvBreakend);
handle!(Pair, LinkedPair);
handle!(Db, DbPair);
handle!(Cluster, SvCluster);
handle!(Chain, SvChain);

/// Whether a list holds this handle (Java's List.contains for classes without equals)
pub fn contains<T: PartialEq>(list: &[T], x: &T) -> bool {
    list.iter().any(|y| y == x)
}
/// List.remove(Object): the first occurrence
pub fn remove_item<T: PartialEq>(list: &mut Vec<T>, x: &T) -> bool {
    if let Some(i) = list.iter().position(|y| y == x) {
        list.remove(i);
        true
    } else {
        false
    }
}
/// List.indexOf
pub fn index_of<T: PartialEq>(list: &[T], x: &T) -> i32 {
    list.iter().position(|y| y == x).map(|i| i as i32).unwrap_or(-1)
}

/// segmentation.Arm
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Arm {
    P,
    Q,
}
impl Arm {
    pub fn as_str(self) -> &'static str {
        match self {
            Arm::P => "P",
            Arm::Q => "Q",
        }
    }
}
impl std::fmt::Display for Arm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

pub fn make_chr_arm_str(chr: &str, arm: Arm) -> String {
    format!("{chr}_{arm}")
}

/// LinxBreakend.coordsStr: chromosome:position:orientation
pub fn coords_str(chromosome: &str, position: i32, orientation: i8) -> String {
    format!("{chromosome}:{position}:{orientation}")
}

/// ResolvedType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ResolvedType {
    NONE,
    DUP_BE,
    LOW_VAF,
    PAIR_INF,
    SGL_MAPPED_INF,
    DEL,
    DUP,
    INF,
    SGL,
    INV,
    INS,
    UNBAL_TRANS,
    SIMPLE_GRP,
    RECIP_TRANS,
    RECIP_TRANS_DEL_DUP,
    RECIP_TRANS_DUPS,
    RECIP_INV,
    RECIP_INV_DUPS,
    RECIP_INV_DEL_DUP,
    DUP_TI,
    DEL_TI,
    UNBAL_TRANS_TI,
    PAIR_OTHER,
    RESOLVED_FOLDBACK,
    FB_INV_PAIR,
    SGL_PAIR_DEL,
    SGL_PAIR_DUP,
    SGL_PAIR_INS,
    LINE,
    COMPLEX,
    DOUBLE_MINUTE,
}

impl ResolvedType {
    pub fn is_simple(self) -> bool {
        use ResolvedType as R;
        matches!(self, R::DEL | R::DUP | R::INS | R::SGL_PAIR_DEL | R::SGL_PAIR_DUP | R::SGL_PAIR_INS)
    }
    pub fn name(self) -> String {
        format!("{self:?}")
    }
}
impl std::fmt::Display for ResolvedType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// ClusteringReason
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ClusteringReason {
    NONE,
    PROXIMITY,
    LOH,
    HOM_LOSS,
    COMMON_ARMS,
    FOLDBACKS,
    CONSEC_BREAKS,
    OVERLAP_FOLDBACKS,
    LOH_CHAIN,
    LONG_DEL_DUP_INV,
    SATELLITE_SGL,
    MAJOR_ALLELE_JCN,
    HIGH_JCN,
    TI_JCN_MATCH,
    SGL_MAPPING_INF,
}

impl ClusteringReason {
    pub const VALUES: [ClusteringReason; 15] = [
        Self::NONE,
        Self::PROXIMITY,
        Self::LOH,
        Self::HOM_LOSS,
        Self::COMMON_ARMS,
        Self::FOLDBACKS,
        Self::CONSEC_BREAKS,
        Self::OVERLAP_FOLDBACKS,
        Self::LOH_CHAIN,
        Self::LONG_DEL_DUP_INV,
        Self::SATELLITE_SGL,
        Self::MAJOR_ALLELE_JCN,
        Self::HIGH_JCN,
        Self::TI_JCN_MATCH,
        Self::SGL_MAPPING_INF,
    ];
    pub fn name(self) -> String {
        format!("{self:?}")
    }
    pub fn value_of(s: &str) -> Option<Self> {
        Self::VALUES.iter().copied().find(|r| r.name() == s)
    }
}
impl std::fmt::Display for ClusteringReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// LineElementType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LineElementType {
    NONE,
    KNOWN,
    SUSPECT,
}

/// A HashSet<LineElementType>: its only reader, LineElementType.toString(set), does not depend on order
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LineElements {
    known: bool,
    suspect: bool,
    none: bool,
}
impl LineElements {
    pub fn add(&mut self, t: LineElementType) {
        match t {
            LineElementType::KNOWN => self.known = true,
            LineElementType::SUSPECT => self.suspect = true,
            LineElementType::NONE => self.none = true,
        }
    }
    pub fn contains(&self, t: LineElementType) -> bool {
        match t {
            LineElementType::KNOWN => self.known,
            LineElementType::SUSPECT => self.suspect,
            LineElementType::NONE => self.none,
        }
    }
    pub fn is_empty(&self) -> bool {
        !self.known && !self.suspect && !self.none
    }
    pub fn len(&self) -> usize {
        self.known as usize + self.suspect as usize + self.none as usize
    }
    /// LineElementType.toString(set)
    pub fn to_str(&self) -> &'static str {
        match self.len() {
            0 => "NONE",
            1 => {
                if self.known {
                    "KNOWN"
                } else if self.suspect {
                    "SUSPECT"
                } else {
                    "NONE"
                }
            }
            _ => "KNOWN;SUSPECT",
        }
    }
}

/// ArmClusterType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ArmClusterType {
    UNKNOWN,
    ISOLATED_BE,
    TI_ONLY,
    DSB,
    FOLDBACK,
    FOLDBACK_DSB,
    COMPLEX_FOLDBACK,
    COMPLEX_LINE,
    COMPLEX_OTHER,
    SIMPLE_DUP,
    SAME_ORIENT,
}
impl std::fmt::Display for ArmClusterType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

// LinxConstants
pub const DEFAULT_PROXIMITY_DISTANCE: i32 = 5000;
pub const MIN_TEMPLATED_INSERTION_LENGTH: i32 = 30;
pub const NO_DB_MARKER: i32 = -1000;
pub const MIN_DEL_LENGTH: i32 = 32;
pub const SHORT_TI_LENGTH: i32 = 1000;
pub const SHORT_DB_LENGTH: i32 = 100;
pub const MAX_MERGE_DISTANCE: i32 = 5000000;
pub const HIGH_JCN_THRESHOLD: f64 = 5.0;
pub const ADJACENT_JCN_RATIO: f64 = 2.3;
pub const MAX_FOLDBACK_CHAIN_LENGTH: i32 = 5000;
pub const MAX_SIMPLE_DUP_DEL_CUTOFF: i32 = 5000000;
pub const MIN_SIMPLE_DUP_DEL_CUTOFF: i32 = 100000;
pub const LOW_JCN_THRESHOLD: f64 = 0.75;
pub const DISRUPTION_CN_REPORTABLE_THRESHOLD: f64 = 0.5;

/// SglMapping: an alternative alignment of a single breakend's insert sequence
#[derive(Clone, Debug)]
pub struct SglMapping {
    pub chromosome: String,
    pub position: i32,
    pub orientation: i8,
    pub cigar: String,
    pub qual_score: i32,
}

impl SglMapping {
    /// SglMapping.from: expected format 4:9973661|-|26S37M11S|19
    pub fn from(mapping_str: &str, breakend_orientation: i8) -> Option<SglMapping> {
        let items: Vec<&str> = mapping_str.split('|').collect();
        if items.len() != 4 {
            return None;
        }
        let location = hmf_common::java_text::split(items[0], ':');
        if location.len() != 2 {
            return None;
        }
        if items[1].is_empty() {
            return None;
        }
        let chromosome = location[0];
        if !hmf_common::sv::is_human_chromosome(chromosome) {
            return None;
        }
        // Java's Integer.parseInt throws on a bad number
        let mut position: i32 = location[1].parse().expect("SGL mapping position");
        let cigar = items[2];
        let strand_sequence: i8 = if items[1] == "+" { 1 } else { -1 };
        let orientation: i8 = if strand_sequence == breakend_orientation { -1 } else { 1 };
        if orientation == 1 {
            position += crate::utils::calc_cigar_aligned_length(cigar) - 1;
        }
        let qual_score: i32 = if !items[3].is_empty() { items[3].parse().expect("SGL mapping qual") } else { 0 };
        Some(SglMapping { chromosome: chromosome.to_string(), position, orientation, cigar: cigar.to_string(), qual_score })
    }

    pub fn convert_from_insert_sequence_alignments(alignments: &str, breakend_orientation: i8) -> Vec<SglMapping> {
        if alignments.is_empty() {
            return Vec::new();
        }
        alignments.split(',').filter_map(|m| SglMapping::from(m, breakend_orientation)).collect()
    }
}

/// SvVarData: one structural variant with everything LINX attaches to it
pub struct SvVarData {
    sv_data: Rc<StructuralVariantData>,
    chr_short: [String; 2],
    arm: Cell<[Arm; 2]>,
    breakend: RefCell<[Option<Be>; 2]>,
    fragile_site: Cell<[bool; 2]>,
    line_elements: Cell<[LineElements; 2]>,
    cluster: RefCell<Option<Cluster>>,
    cluster_reason: RefCell<String>,
    foldback_breakends: RefCell<[Option<Be>; 2]>,
    foldback_length: Cell<[i32; 2]>,
    foldback_info: RefCell<[String; 2]>,
    nearest_sv_distance: Cell<i32>,
    nearest_sv_relation: RefCell<String>,
    ti_links: RefCell<[Vec<Pair>; 2]>,
    db_link: RefCell<[Option<Db>; 2]>,
    assembly_links: RefCell<[Vec<String>; 2]>,
    genes: RefCell<[Vec<BreakendGene>; 2]>,
    copy_number: Cell<[f64; 2]>,
    copy_number_change: Cell<[f64; 2]>,
    jcn: Cell<f64>,
    has_calc_jcn: Cell<bool>,
    jcn_min: Cell<f64>,
    jcn_max: Cell<f64>,
    /// prev start, post start, prev end, post end
    cn_data: RefCell<[Option<SvCnData>; 4]>,
    sgl_mappings: RefCell<Vec<SglMapping>>,
    linked_svs: RefCell<Option<[Sv; 2]>>,
}

pub const RELATION_TYPE_NEIGHBOUR: &str = "NHBR";
pub const RELATION_TYPE_OVERLAP: &str = "OVRL";
pub const SGL_CENTRO_SATELLITE: &str = "Satellite/centr";
pub const SGL_TELO_SATELLITE: &str = "Satellite/telo";
pub const INF_SV_TYPE: &str = "INF";
pub const CR_DELIM: &str = "-";
pub const PRE_TRANSCRIPT_DISTANCE: i32 = 10000;

impl SvVarData {
    pub fn new(sv_data: Rc<StructuralVariantData>) -> SvVarData {
        let strip = |c: &str| hmf_common::regions::strip_chr_prefix(c).to_string();
        let chr_short = [strip(&sv_data.start_chromosome), strip(&sv_data.end_chromosome)];
        let copy_number = [sv_data.adjusted_start_copy_number, sv_data.adjusted_end_copy_number];
        let copy_number_change = [sv_data.adjusted_start_copy_number_change, sv_data.adjusted_end_copy_number_change];
        let jcn = sv_data.junction_copy_number;
        let assembly = [assembly_data(&sv_data.start_linked_by), assembly_data(&sv_data.end_linked_by)];
        let t = sv_data.sv_type();
        let sgl_mappings = if t == SGL || t == INF {
            SglMapping::convert_from_insert_sequence_alignments(&sv_data.insert_sequence_alignments, sv_data.start_orientation)
        } else {
            Vec::new()
        };
        SvVarData {
            sv_data,
            chr_short,
            arm: Cell::new([Arm::P, Arm::P]),
            breakend: RefCell::new([None, None]),
            fragile_site: Cell::new([false, false]),
            line_elements: Cell::new([LineElements::default(), LineElements::default()]),
            cluster: RefCell::new(None),
            cluster_reason: RefCell::new(String::new()),
            foldback_breakends: RefCell::new([None, None]),
            foldback_length: Cell::new([-1, -1]),
            foldback_info: RefCell::new([String::new(), String::new()]),
            nearest_sv_distance: Cell::new(-1),
            nearest_sv_relation: RefCell::new(String::new()),
            ti_links: RefCell::new([Vec::new(), Vec::new()]),
            db_link: RefCell::new([None, None]),
            assembly_links: RefCell::new(assembly),
            genes: RefCell::new([Vec::new(), Vec::new()]),
            copy_number: Cell::new(copy_number),
            copy_number_change: Cell::new(copy_number_change),
            jcn: Cell::new(jcn),
            has_calc_jcn: Cell::new(false),
            jcn_min: Cell::new(0.0),
            jcn_max: Cell::new(0.0),
            cn_data: RefCell::new([None, None, None, None]),
            sgl_mappings: RefCell::new(sgl_mappings),
            linked_svs: RefCell::new(None),
        }
    }

    pub fn id(&self) -> i32 {
        self.sv_data.id
    }
    pub fn sv_data(&self) -> &StructuralVariantData {
        &self.sv_data
    }
    pub fn sv_data_rc(&self) -> &Rc<StructuralVariantData> {
        &self.sv_data
    }
    pub fn chromosome(&self, is_start: bool) -> &str {
        if is_start { &self.sv_data.start_chromosome } else { &self.sv_data.end_chromosome }
    }
    pub fn position(&self, is_start: bool) -> i32 {
        if is_start { self.sv_data.start_position } else { self.sv_data.end_position }
    }
    pub fn orientation(&self, is_start: bool) -> i8 {
        if is_start { self.sv_data.start_orientation } else { self.sv_data.end_orientation }
    }
    pub fn copy_number(&self, is_start: bool) -> f64 {
        self.copy_number.get()[se(is_start)]
    }
    pub fn sv_type(&self) -> StructuralVariantType {
        self.sv_data.sv_type()
    }
    pub fn breakend(&self, is_start: bool) -> Option<Be> {
        self.breakend.borrow()[se(is_start)].clone()
    }
    /// getBreakend for a breakend known to exist
    pub fn be(&self, is_start: bool) -> Be {
        self.breakend(is_start).expect("SV breakend")
    }
    pub fn breakend_se(&self, se: usize) -> Option<Be> {
        self.breakend.borrow()[se].clone()
    }
    pub fn set_sgl_end_breakend(&self, breakend: Option<Be>) {
        self.breakend.borrow_mut()[SE_END] = breakend;
    }
    pub fn is_sgl_breakend(&self) -> bool {
        let t = self.sv_type();
        t == SGL || t == INF
    }
    pub fn is_inferred_sgl(&self) -> bool {
        self.sv_type() == INF
    }

    pub fn pos_id(&self) -> String {
        if self.is_sgl_breakend() {
            format!("id({}) pos({})", self.id(), coords_str(&self.chr_short[0], self.position(true), self.orientation(true)))
        } else {
            format!(
                "id({}) pos({} -> {})",
                self.id(),
                coords_str(&self.chr_short[0], self.position(true), self.orientation(true)),
                coords_str(&self.chr_short[1], self.position(false), self.orientation(false))
            )
        }
    }

    pub fn coords_str(&self, use_start: bool) -> String {
        if self.is_sgl_breakend() && !use_start {
            return String::new();
        }
        coords_str(self.chromosome(use_start), self.position(use_start), self.orientation(use_start))
    }

    pub fn arm(&self, is_start: bool) -> Arm {
        self.arm.get()[se(is_start)]
    }
    pub fn chr_short(&self, is_start: bool) -> &str {
        &self.chr_short[se(is_start)]
    }

    pub fn cluster(&self) -> Option<Cluster> {
        self.cluster.borrow().clone()
    }
    /// getCluster for an SV known to be clustered
    pub fn cl(&self) -> Cluster {
        self.cluster().expect("SV cluster")
    }
    pub fn set_cluster(&self, cluster: Option<Cluster>) {
        *self.cluster.borrow_mut() = cluster;
    }

    pub fn length(&self) -> i32 {
        if self.sv_type() == BND || self.is_sgl_breakend() {
            return 0;
        }
        (self.position(false) - self.position(true)).abs()
    }

    pub fn cluster_reason(&self) -> String {
        self.cluster_reason.borrow().clone()
    }
    pub fn has_cluster_reason(&self, reason: ClusteringReason) -> bool {
        self.cluster_reason.borrow().contains(&reason.name())
    }

    pub fn jcn(&self) -> f64 {
        self.jcn.get()
    }

    pub fn implied_jcn(&self) -> i32 {
        self.max_assembled_breakend().max(self.rounded_jcn(true) as i32)
    }

    pub fn copy_number_change(&self, is_start: bool) -> f64 {
        if let Some(db) = &self.db_link.borrow()[se(is_start)] {
            if db.length() == 0 {
                return self.jcn.get();
            }
        }
        self.copy_number_change.get()[se(is_start)]
    }

    pub fn set_copy_number_data(&self, is_start: bool, copy_number: f64, copy_number_change: f64) {
        let mut cn = self.copy_number.get();
        cn[se(is_start)] = copy_number;
        self.copy_number.set(cn);
        let mut cc = self.copy_number_change.get();
        cc[se(is_start)] = copy_number_change;
        self.copy_number_change.set(cc);
    }

    pub fn rounded_jcn(&self, enforce_clonal: bool) -> f64 {
        let r = crate::utils::jround(self.jcn()) as f64;
        if enforce_clonal { r.max(1.0) } else { r }
    }

    pub fn max_assembled_breakend(&self) -> i32 {
        let l = self.ti_links.borrow();
        let a = l[0].iter().filter(|p| p.is_assembled()).count();
        let b = l[1].iter().filter(|p| p.is_assembled()).count();
        a.max(b) as i32
    }

    pub fn nearest_sv_distance(&self) -> i32 {
        self.nearest_sv_distance.get()
    }
    pub fn set_nearest_sv_distance(&self, d: i32) {
        self.nearest_sv_distance.set(d)
    }
    pub fn nearest_sv_relation(&self) -> String {
        self.nearest_sv_relation.borrow().clone()
    }
    pub fn set_nearest_sv_relation(&self, r: &str) {
        *self.nearest_sv_relation.borrow_mut() = r.to_string();
    }

    pub fn set_fragile_sites(&self, start: bool, end: bool) {
        self.fragile_site.set([start, end]);
    }
    pub fn is_fragile_site(&self, use_start: bool) -> bool {
        self.fragile_site.get()[se(use_start)]
    }

    pub fn add_line_element(&self, t: LineElementType, is_start: bool) {
        let mut l = self.line_elements.get();
        l[se(is_start)].add(t);
        self.line_elements.set(l);
    }
    pub fn is_line_element(&self, use_start: bool) -> bool {
        !self.line_elements.get()[se(use_start)].is_empty()
    }
    pub fn in_line_element(&self) -> bool {
        self.is_line_element(true) || self.is_line_element(false)
    }
    pub fn has_line_element(&self, t: LineElementType, is_start: bool) -> bool {
        self.line_elements.get()[se(is_start)].contains(t)
    }
    pub fn line_element(&self, use_start: bool) -> LineElements {
        self.line_elements.get()[se(use_start)]
    }

    pub fn linked_pairs(&self, is_start: bool) -> Vec<Pair> {
        self.ti_links.borrow()[se(is_start)].clone()
    }
    pub fn assembled_linked_pairs(&self, is_start: bool) -> Vec<Pair> {
        self.ti_links.borrow()[se(is_start)].iter().filter(|p| p.is_assembled()).cloned().collect()
    }
    pub fn linked_pair(&self, is_start: bool) -> Option<Pair> {
        self.ti_links.borrow()[se(is_start)].first().cloned()
    }

    /// addLinkedPair: in order from shortest to longest, without duplicates
    pub fn add_linked_pair(&self, link: &Pair, is_start: bool) {
        let mut links = self.ti_links.borrow_mut();
        let links = &mut links[se(is_start)];
        let mut index = 0;
        while index < links.len() {
            let other = &links[index];
            if other.matches(link) {
                return;
            }
            if link.position_distance() < links[index].position_distance() {
                break;
            }
            index += 1;
        }
        links.insert(index, link.clone());
    }

    /// getLinkedPairs(isStart).remove(pair)
    pub fn remove_linked_pair(&self, pair: &Pair, is_start: bool) {
        remove_item(&mut self.ti_links.borrow_mut()[se(is_start)], pair);
    }

    pub fn clear_clustering_data(&self) {
        self.cluster_reason.borrow_mut().clear();
        let mut l = self.ti_links.borrow_mut();
        l[0].clear();
        l[1].clear();
    }

    pub fn db_link(&self, is_start: bool) -> Option<Db> {
        self.db_link.borrow()[se(is_start)].clone()
    }
    pub fn set_db_link(&self, link: Option<Db>, is_start: bool) {
        self.db_link.borrow_mut()[se(is_start)] = link;
    }

    pub fn foldback_id(&self, use_start: bool) -> i32 {
        match &self.foldback_breakends.borrow()[se(use_start)] {
            Some(b) => b.sv().id(),
            None => -1,
        }
    }
    pub fn foldback_breakend(&self, is_start: bool) -> Option<Be> {
        self.foldback_breakends.borrow()[se(is_start)].clone()
    }
    pub fn foldback_length(&self, is_start: bool) -> i32 {
        self.foldback_length.get()[se(is_start)]
    }
    pub fn foldback_info(&self, is_start: bool) -> String {
        self.foldback_info.borrow()[se(is_start)].clone()
    }
    pub fn is_foldback(&self) -> bool {
        let f = self.foldback_breakends.borrow();
        f[0].is_some() || f[1].is_some()
    }

    pub fn is_single_breakend_foldback(&self) -> bool {
        let f = self.foldback_breakends.borrow();
        let b = self.breakend.borrow();
        (f[0].is_some() && f[0] == b[0]) || (f[1].is_some() && f[1] == b[1])
    }

    pub fn is_chained_foldback(&self) -> bool {
        let f = self.foldback_breakends.borrow();
        let b = self.breakend.borrow();
        (f[1].is_some() && f[1] != b[0]) || (f[0].is_some() && f[0] != b[1])
    }

    pub fn type_str(&self) -> &'static str {
        if self.is_inferred_sgl() { INF_SV_TYPE } else { self.sv_type().name() }
    }

    /// both ends within the same chromosomal arm
    pub fn is_local(&self) -> bool {
        self.chromosome(true) == self.chromosome(false) && self.arm(true) == self.arm(false)
    }
    pub fn is_cross_arm(&self) -> bool {
        !self.is_sgl_breakend() && !self.is_local()
    }
    pub fn is_simple_type(&self) -> bool {
        matches!(self.sv_type(), DEL | DUP | INS)
    }

    pub fn linked_svs(&self) -> Option<[Sv; 2]> {
        self.linked_svs.borrow().clone()
    }
    pub fn set_linked_svs(&self, var1: &Sv, var2: &Sv) {
        *self.linked_svs.borrow_mut() = Some([var1.clone(), var2.clone()]);
    }

    pub fn genes_list(&self, is_start: bool) -> Vec<BreakendGene> {
        self.genes.borrow()[se(is_start)].clone()
    }
    pub fn genes_mut(&self) -> RefMut<'_, [Vec<BreakendGene>; 2]> {
        self.genes.borrow_mut()
    }
    pub fn genes_ref(&self) -> Ref<'_, [Vec<BreakendGene>; 2]> {
        self.genes.borrow()
    }

    pub fn sgl_mappings(&self) -> Vec<SglMapping> {
        self.sgl_mappings.borrow().clone()
    }
    pub fn has_sgl_mappings(&self) -> bool {
        !self.sgl_mappings.borrow().is_empty()
    }
    pub fn add_sgl_mapping(&self, m: SglMapping) {
        self.sgl_mappings.borrow_mut().push(m);
    }

    /// getGeneInBreakend(isStart, includeId, includeTransImpact)
    pub fn gene_in_breakend(&self, is_start: bool, include_id: bool, include_trans_impact: bool) -> String {
        let genes = self.genes.borrow();
        let mut items: Vec<String> = Vec::new();
        for gene in genes[se(is_start)].iter().filter(|g| g.breakend_within_gene(PRE_TRANSCRIPT_DISTANCE)) {
            let mut g: Vec<String> = Vec::new();
            if include_id {
                g.push(gene.gene_id().to_string());
            }
            g.push(gene.gene_name().to_string());
            if include_trans_impact {
                if let Some(trans) = gene.canonical() {
                    g.push(trans.region_type().to_string());
                    g.push(trans.coding_type().to_string());
                    if trans.is_disruptive() {
                        g.push("disruptive".to_string());
                    }
                    if trans.coding_type() == crate::gene::TranscriptCodingType::CODING {
                        g.push(format!("exon={}", trans.exon_upstream));
                    }
                }
            }
            items.push(g.join("|"));
        }
        items.join(";")
    }

    pub fn has_assembly_link(&self, is_start: bool) -> bool {
        self.ti_links.borrow()[se(is_start)].iter().any(|p| p.is_assembled())
    }
    pub fn assembly_links(&self, is_start: bool) -> Vec<String> {
        self.assembly_links.borrow()[se(is_start)].clone()
    }
    pub fn assembly_info_str(&self, is_start: bool) -> String {
        self.assembly_links.borrow()[se(is_start)].join(";")
    }
    pub fn add_assembly_info(&self, is_start: bool, data: &str) {
        self.assembly_links.borrow_mut()[se(is_start)].push(data.to_string());
    }

    pub fn set_jcn_recalc_data(&self, min_jcn: f64, max_jcn: f64) {
        if max_jcn == 0.0 {
            // suggests an error in the calculation
            self.jcn_max.set(self.jcn.get());
        } else {
            self.jcn_min.set(min_jcn);
            self.jcn_max.set(max_jcn);
            self.jcn.set((min_jcn + max_jcn) * 0.5);
        }
        self.has_calc_jcn.set(true);
    }
    pub fn has_calculated_jcn(&self) -> bool {
        self.has_calc_jcn.get()
    }
    pub fn jcn_max(&self) -> f64 {
        if self.has_calc_jcn.get() { self.jcn_max.get() } else { self.jcn.get() }
    }
    pub fn jcn_min(&self) -> f64 {
        if self.has_calc_jcn.get() { self.jcn_min.get() } else { self.jcn.get() }
    }
    pub fn jcn_uncertainty(&self) -> f64 {
        if self.has_calc_jcn.get() { (self.jcn.get() - self.jcn_min.get()) * 0.5 } else { 0.0 }
    }

    pub fn calc_vaf(&self, is_start: bool) -> f64 {
        self.copy_number_change(is_start) / self.copy_number(is_start)
    }

    pub fn copy_number_data(&self, is_start: bool, is_previous: bool) -> Option<SvCnData> {
        let i = if is_start { if is_previous { 0 } else { 1 } } else if is_previous { 2 } else { 3 };
        self.cn_data.borrow()[i].clone()
    }
    pub fn set_copy_number_segments(&self, is_start: bool, prev: Option<SvCnData>, post: Option<SvCnData>) {
        let mut c = self.cn_data.borrow_mut();
        if is_start {
            c[0] = prev;
            c[1] = post;
        } else {
            c[2] = prev;
            c[3] = post;
        }
    }

    pub fn sgl_to_centromere_or_telomere(&self) -> bool {
        if !self.is_sgl_breakend() || self.is_inferred_sgl() {
            return false;
        }
        let rc = &self.sv_data.insert_sequence_repeat_class;
        if rc == SGL_CENTRO_SATELLITE || rc == SGL_TELO_SATELLITE {
            return true;
        }
        if rc == "Simple_repeat" {
            let rt = &self.sv_data.insert_sequence_repeat_type;
            if rt == "(CCCTAA)n" || rt == "(TTAGGG)n" {
                return true;
            }
        }
        false
    }

    pub fn sgl_to_satellite_repeats(&self) -> bool {
        if self.sv_type() != SGL {
            return false;
        }
        if self.sv_data.insert_sequence_repeat_class == SGL_CENTRO_SATELLITE {
            return true;
        }
        let rt = self.sv_data.insert_sequence_repeat_type.as_str();
        rt == "(CATTC)n" || rt == "(GAATG)n" || rt == "HSATII" || rt == "SAR"
    }
}

fn assembly_data(data: &str) -> Vec<String> {
    if data.is_empty() || data == "." {
        return Vec::new();
    }
    let delim = if data.contains(';') { ';' } else { ',' };
    data.split(delim).map(|s| s.to_string()).collect()
}

impl Sv {
    pub fn from_data(data: Rc<StructuralVariantData>) -> Sv {
        Sv::new(SvVarData::new(data))
    }

    /// setChromosomalArms: also creates the breakends
    pub fn set_chromosomal_arms(&self, start: Arm, end: Arm) {
        self.arm.set([start, end]);
        let s = Be::new(SvBreakend::new(self, true));
        let e = if !self.is_sgl_breakend() { Some(Be::new(SvBreakend::new(self, false))) } else { None };
        let mut b = self.breakend.borrow_mut();
        b[SE_START] = Some(s);
        if let Some(e) = e {
            b[SE_END] = Some(e);
        }
    }

    /// addClusterReason
    pub fn add_cluster_reason(&self, reason: ClusteringReason, other_id: i32) {
        let name = reason.name();
        {
            let mut r = self.cluster_reason.borrow_mut();
            if r.contains(&name) {
                return;
            }
            crate::utils::append_str(&mut r, &name, ';');
            if other_id > -1 {
                r.push_str(CR_DELIM);
                r.push_str(&other_id.to_string());
            }
        }
        if other_id == self.id() {
            hmf_common::warn!("SV({}) reason({}) setting to own ID", self.id(), reason);
        }
    }

    /// setFoldbackLink
    pub fn set_foldback_link(&self, is_start: bool, link: Option<Be>, length: i32, link_info: &str) {
        self.foldback_breakends.borrow_mut()[se(is_start)] = link;
        let mut l = self.foldback_length.get();
        l[se(is_start)] = length;
        self.foldback_length.set(l);
        self.foldback_info.borrow_mut()[se(is_start)] = link_info.to_string();
        if let Some(cluster) = self.cluster() {
            let none = {
                let f = self.foldback_breakends.borrow();
                f[0].is_none() && f[1].is_none()
            };
            if none {
                cluster.deregister_foldback(self);
            } else {
                cluster.register_foldback(self);
            }
        }
    }

    pub fn break_links(&self) {
        *self.cluster.borrow_mut() = None;
        *self.breakend.borrow_mut() = [None, None];
        *self.foldback_breakends.borrow_mut() = [None, None];
        *self.ti_links.borrow_mut() = [Vec::new(), Vec::new()];
        *self.db_link.borrow_mut() = [None, None];
        *self.linked_svs.borrow_mut() = None;
        *self.genes.borrow_mut() = [Vec::new(), Vec::new()];
        *self.cn_data.borrow_mut() = [None, None, None, None];
    }
}

impl std::fmt::Display for Sv {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.pos_id(), self.type_str())
    }
}

/// SvVarData.haveSameChrArms
pub fn have_same_chr_arms(var1: &SvVarData, var2: &SvVarData) -> bool {
    (var1.chromosome(true) == var2.chromosome(true)
        && var1.chromosome(false) == var2.chromosome(false)
        && var1.arm(true) == var2.arm(true)
        && var1.arm(false) == var2.arm(false))
        || (var1.chromosome(true) == var2.chromosome(false)
            && var1.chromosome(false) == var2.chromosome(true)
            && var1.arm(true) == var2.arm(false)
            && var1.arm(false) == var2.arm(true))
}

/// SvBreakend
pub struct SvBreakend {
    sv: Sv,
    chr_pos_index: Cell<i32>,
    cluster_chr_pos_index: Cell<i32>,
    chromosome: String,
    arm: Arm,
    chr_arm: String,
    position: i32,
    orientation: i8,
    uses_start: bool,
}

pub const DIRECTION_CENTROMERE: &str = "C";
pub const DIRECTION_TELOMERE: &str = "T";

impl SvBreakend {
    pub fn new(var: &Sv, use_start: bool) -> SvBreakend {
        let chromosome = var.chromosome(use_start).to_string();
        let arm = var.arm(use_start);
        SvBreakend {
            sv: var.clone(),
            chr_pos_index: Cell::new(-1),
            cluster_chr_pos_index: Cell::new(-1),
            chr_arm: make_chr_arm_str(&chromosome, arm),
            chromosome,
            arm,
            position: var.position(use_start),
            orientation: var.orientation(use_start),
            uses_start: use_start,
        }
    }

    /// a breakend at an SGL's mapping location
    pub fn from_mapping(var: &Sv, mapping: &SglMapping) -> SvBreakend {
        let arm = crate::utils::get_chromosomal_arm(&mapping.chromosome, mapping.position);
        SvBreakend {
            sv: var.clone(),
            chr_pos_index: Cell::new(-1),
            cluster_chr_pos_index: Cell::new(-1),
            chr_arm: make_chr_arm_str(&mapping.chromosome, arm),
            chromosome: mapping.chromosome.clone(),
            arm,
            position: mapping.position,
            orientation: mapping.orientation,
            uses_start: false,
        }
    }

    pub fn sv(&self) -> &Sv {
        &self.sv
    }
    pub fn cluster(&self) -> Option<Cluster> {
        self.sv.cluster()
    }
    pub fn cl(&self) -> Cluster {
        self.sv.cl()
    }
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }
    pub fn arm(&self) -> Arm {
        self.arm
    }
    pub fn chr_arm(&self) -> &str {
        &self.chr_arm
    }
    pub fn position(&self) -> i32 {
        self.position
    }
    pub fn orientation(&self) -> i8 {
        self.orientation
    }
    pub fn forward_orient(&self) -> bool {
        self.orientation == 1
    }
    pub fn reverse_orient(&self) -> bool {
        self.orientation == -1
    }
    pub fn uses_start(&self) -> bool {
        self.uses_start
    }
    pub fn sv_type(&self) -> StructuralVariantType {
        self.sv.sv_type()
    }
    pub fn other_breakend(&self) -> Option<Be> {
        self.sv.breakend(!self.uses_start)
    }
    pub fn set_chr_pos_index(&self, i: i32) {
        self.chr_pos_index.set(i)
    }
    pub fn chr_pos_index(&self) -> i32 {
        self.chr_pos_index.get()
    }
    pub fn set_cluster_chr_pos_index(&self, i: i32) {
        self.cluster_chr_pos_index.set(i)
    }
    pub fn cluster_chr_pos_index(&self) -> i32 {
        self.cluster_chr_pos_index.get()
    }
    pub fn direction(&self) -> &'static str {
        if self.forward_orient() == (self.arm == Arm::P) { DIRECTION_TELOMERE } else { DIRECTION_CENTROMERE }
    }
    pub fn anchor_distance(&self) -> i32 {
        if !self.uses_start && self.sv.is_sgl_breakend() {
            return 0;
        }
        if self.uses_start {
            self.sv.sv_data().start_anchoring_support_distance
        } else {
            self.sv.sv_data().end_anchoring_support_distance
        }
    }
    pub fn copy_number_low_side(&self) -> f64 {
        self.sv.copy_number(self.uses_start) - self.sv.copy_number_change(self.uses_start)
    }
    pub fn copy_number_change(&self) -> f64 {
        self.sv.copy_number_change(self.uses_start)
    }
    pub fn copy_number(&self) -> f64 {
        self.sv.copy_number(self.uses_start)
    }
    pub fn jcn(&self) -> f64 {
        self.sv.jcn()
    }
    pub fn jcn_uncertainty(&self) -> f64 {
        self.sv.jcn_uncertainty()
    }
    pub fn homology(&self) -> String {
        if !self.uses_start && self.sv.is_sgl_breakend() {
            return String::new();
        }
        if self.uses_start {
            self.sv.sv_data().start_homology_sequence.clone()
        } else {
            self.sv.sv_data().end_homology_sequence.clone()
        }
    }
    /// getCopyNumber(usePrevious)
    pub fn get_copy_number(&self, use_previous: bool) -> f64 {
        self.sv.copy_number_data(self.uses_start, use_previous).map(|c| c.copy_number).unwrap_or(0.0)
    }
    pub fn major_allele_jcn(&self, use_previous: bool) -> f64 {
        self.sv.copy_number_data(self.uses_start, use_previous).map(|c| c.major_allele_jcn()).unwrap_or(0.0)
    }
    pub fn minor_allele_jcn(&self, use_previous: bool) -> f64 {
        self.sv.copy_number_data(self.uses_start, use_previous).map(|c| c.minor_allele_jcn()).unwrap_or(0.0)
    }
    pub fn is_assembled_link(&self) -> bool {
        self.sv.has_assembly_link(self.uses_start)
    }
    pub fn db_link(&self) -> Option<Db> {
        self.sv.db_link(self.uses_start)
    }
    pub fn linked_pairs(&self) -> Vec<Pair> {
        self.sv.linked_pairs(self.uses_start)
    }
    pub fn is_foldback(&self) -> bool {
        self.sv.foldback_breakend(self.uses_start).is_some()
    }
    pub fn foldback_breakend(&self) -> Option<Be> {
        self.sv.foldback_breakend(self.uses_start)
    }
    pub fn in_line_element(&self) -> bool {
        self.sv.is_line_element(self.uses_start)
    }
    pub fn has_line_element(&self, t: LineElementType) -> bool {
        self.sv.has_line_element(t, self.uses_start)
    }
    pub fn genes_list(&self) -> Vec<BreakendGene> {
        self.sv.genes_list(self.uses_start)
    }
}

impl std::fmt::Display for SvBreakend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: {} {}:{}:{}",
            self.sv.id(),
            if self.uses_start { "start" } else { "end" },
            self.chromosome,
            self.position,
            self.orientation
        )
    }
}
impl std::fmt::Display for Be {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// LinkedPair: a templated insertion between two breakends
pub struct LinkedPair {
    first: RefCell<Be>,
    second: RefCell<Be>,
    link_length: i32,
    is_inferred: Cell<bool>,
    link_reason: RefCell<String>,
    link_index: Cell<i32>,
    traversed_sv_count: Cell<i32>,
    next_sv_distance: Cell<i32>,
    next_clustered_sv_distance: Cell<i32>,
    location_type: RefCell<&'static str>,
    overlap_count: Cell<i32>,
    has_copy_number_gain: Cell<bool>,
    exon_match_data: RefCell<String>,
}

pub const LOCATION_TYPE_UNCLEAR: &str = "Unclear";
pub const LOCATION_TYPE_REMOTE: &str = "Remote";
pub const LOCATION_TYPE_EXTERNAL: &str = "External";
pub const LOCATION_TYPE_INTERNAL: &str = "Internal";

impl LinkedPair {
    pub fn first_breakend(&self) -> Be {
        self.first.borrow().clone()
    }
    pub fn second_breakend(&self) -> Be {
        self.second.borrow().clone()
    }
    pub fn first(&self) -> Sv {
        self.first.borrow().sv().clone()
    }
    pub fn second(&self) -> Sv {
        self.second.borrow().sv().clone()
    }
    pub fn first_link_on_start(&self) -> bool {
        self.first.borrow().uses_start()
    }
    pub fn second_link_on_start(&self) -> bool {
        self.second.borrow().uses_start()
    }
    pub fn first_unlinked_on_start(&self) -> bool {
        !self.first.borrow().uses_start()
    }
    pub fn second_unlinked_on_start(&self) -> bool {
        !self.second.borrow().uses_start()
    }

    /// getBreakend(isStart): the lower (start) or upper breakend; an SGL linked on its start gives its
    /// only breakend
    pub fn breakend(&self, is_start: bool) -> Be {
        let f = self.first_breakend();
        let s = self.second_breakend();
        let be_first = if f.sv().is_sgl_breakend() && f.uses_start() { f.sv().be(true) } else { f };
        let be_second = if s.sv().is_sgl_breakend() && s.uses_start() { s.sv().be(true) } else { s };
        if is_start {
            if be_first.position() <= be_second.position() { be_first } else { be_second }
        } else if be_first.position() > be_second.position() {
            be_first
        } else {
            be_second
        }
    }
    pub fn breakend_se(&self, se: usize) -> Be {
        self.breakend(is_start(se))
    }

    pub fn chromosome(&self) -> String {
        self.first.borrow().chromosome().to_string()
    }
    pub fn position_distance(&self) -> i32 {
        self.link_length
    }
    pub fn base_length(&self) -> i32 {
        self.link_length + 1
    }
    pub fn set_is_assembled(&self) {
        self.is_inferred.set(false)
    }
    pub fn is_inferred(&self) -> bool {
        self.is_inferred.get()
    }
    pub fn is_assembled(&self) -> bool {
        !self.is_inferred.get()
    }
    pub fn link_reason(&self) -> String {
        self.link_reason.borrow().clone()
    }
    pub fn set_link_reason(&self, reason: &str, index: i32) {
        *self.link_reason.borrow_mut() = reason.to_string();
        self.link_index.set(index);
    }
    pub fn link_index(&self) -> i32 {
        self.link_index.get()
    }
    pub fn set_next_sv_data(&self, distance: i32, cluster_distance: i32) {
        self.next_sv_distance.set(distance);
        self.next_clustered_sv_distance.set(cluster_distance);
    }
    pub fn next_sv_distance(&self) -> i32 {
        self.next_sv_distance.get()
    }
    pub fn next_clustered_sv_distance(&self) -> i32 {
        self.next_clustered_sv_distance.get()
    }
    pub fn set_location_type(&self, t: &'static str) {
        *self.location_type.borrow_mut() = t;
    }
    pub fn location_type(&self) -> &'static str {
        *self.location_type.borrow()
    }
    pub fn set_overlap_count(&self, c: i32) {
        self.overlap_count.set(c)
    }
    pub fn overlap_count(&self) -> i32 {
        self.overlap_count.get()
    }
    pub fn set_has_copy_number_gain(&self, t: bool) {
        self.has_copy_number_gain.set(t)
    }
    pub fn has_copy_number_gain(&self) -> bool {
        self.has_copy_number_gain.get()
    }
    pub fn set_traversed_sv_count(&self, c: i32) {
        self.traversed_sv_count.set(c)
    }
    pub fn traversed_sv_count(&self) -> i32 {
        self.traversed_sv_count.get()
    }
    pub fn set_exon_match_data(&self, d: &str) {
        *self.exon_match_data.borrow_mut() = d.to_string();
    }
    pub fn exon_match_data(&self) -> String {
        self.exon_match_data.borrow().clone()
    }

    pub fn has_breakend_of(&self, var: &Sv, use_start: bool) -> bool {
        let f = self.first.borrow();
        let s = self.second.borrow();
        (f.sv() == var && f.uses_start() == use_start) || (s.sv() == var && s.uses_start() == use_start)
    }
    pub fn has_breakend(&self, breakend: &Be) -> bool {
        self.has_breakend_of(breakend.sv(), breakend.uses_start())
    }
    pub fn has_link_clash(&self, other: &LinkedPair) -> bool {
        self.has_breakend_of(&other.first(), other.first_link_on_start()) || self.has_breakend_of(&other.second(), other.second_link_on_start())
    }

    pub fn switch_svs(&self) {
        std::mem::swap(&mut *self.first.borrow_mut(), &mut *self.second.borrow_mut());
    }

    /// matches: the same two breakends, in either order
    pub fn matches(&self, other: &LinkedPair) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        }
        let (f, s) = (self.first_breakend(), self.second_breakend());
        let (of, os) = (other.first_breakend(), other.second_breakend());
        (f == of && s == os) || (f == os && s == of)
    }

    pub fn opposite_match(&self, other: &LinkedPair) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        }
        let (f, s) = (self.first_breakend(), self.second_breakend());
        if f.sv() == &other.first()
            && f.uses_start() != other.first_link_on_start()
            && s.sv() == &other.second()
            && s.uses_start() != other.second_link_on_start()
        {
            return true;
        }
        f.sv() == &other.second()
            && f.uses_start() != other.second_link_on_start()
            && s.sv() == &other.first()
            && s.uses_start() != other.first_link_on_start()
    }

    pub fn same_variants(&self, other: &LinkedPair) -> bool {
        let (a, b) = (self.first(), self.second());
        let (c, d) = (other.first(), other.second());
        (a == c && b == d) || (a == d && b == c)
    }

    pub fn other_sv(&self, var: &Sv) -> Sv {
        if &self.first() == var { self.second() } else { self.first() }
    }
    pub fn other_breakend(&self, breakend: &Be) -> Be {
        let lower = self.breakend(true);
        if breakend == &lower { self.breakend(false) } else { lower }
    }
    pub fn has_variant(&self, var: &Sv) -> bool {
        &self.first() == var || &self.second() == var
    }
    pub fn is_dup_link(&self) -> bool {
        self.first() == self.second() && self.first_breakend().sv_type() == DUP
    }
}

impl Pair {
    pub fn from_breakends(first: &Be, second: &Be) -> Pair {
        let length = (first.position() - second.position()).abs();
        Pair::new(LinkedPair {
            first: RefCell::new(first.clone()),
            second: RefCell::new(second.clone()),
            link_length: length,
            is_inferred: Cell::new(true),
            link_reason: RefCell::new(String::new()),
            link_index: Cell::new(-1),
            traversed_sv_count: Cell::new(0),
            next_sv_distance: Cell::new(0),
            next_clustered_sv_distance: Cell::new(0),
            location_type: RefCell::new(LOCATION_TYPE_UNCLEAR),
            overlap_count: Cell::new(0),
            has_copy_number_gain: Cell::new(false),
            exon_match_data: RefCell::new(String::new()),
        })
    }
    pub fn from_svs(first: &Sv, second: &Sv, first_link_on_start: bool, second_link_on_start: bool) -> Pair {
        Pair::from_breakends(&first.be(first_link_on_start), &second.be(second_link_on_start))
    }
    /// LinkedPair.copy
    pub fn copy_of(other: &Pair) -> Pair {
        let p = Pair::from_breakends(&other.first_breakend(), &other.second_breakend());
        p.set_link_reason(&other.link_reason(), other.link_index());
        if other.is_assembled() {
            p.set_is_assembled();
        }
        p
    }
}

impl std::fmt::Display for LinkedPair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (a, b) = (self.first_breakend(), self.second_breakend());
        write!(
            f,
            "{} {}:{}:{} & {} {}:{}:{}",
            a.sv().id(),
            a.chromosome(),
            a.position(),
            if a.uses_start() { "start" } else { "end" },
            b.sv().id(),
            b.chromosome(),
            b.position(),
            if b.uses_start() { "start" } else { "end" }
        )
    }
}
impl std::fmt::Display for Pair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// LinkedPair.hasLinkClash(links1, links2)
pub fn links_have_clash(links1: &[Pair], links2: &[Pair]) -> bool {
    links1.iter().any(|p| links2.iter().any(|x| p.has_link_clash(x)))
}

/// DbPair: a deletion bridge between two breakends
pub struct DbPair {
    lower: Be,
    upper: Be,
    link_length: i32,
}

impl DbPair {
    pub fn new(breakend1: &Be, breakend2: &Be) -> DbPair {
        let (lower, upper) = if breakend1.position() < breakend2.position()
            || (breakend1.position() == breakend2.position() && breakend1.forward_orient())
        {
            (breakend1.clone(), breakend2.clone())
        } else {
            (breakend2.clone(), breakend1.clone())
        };
        // exact base is considered a 1 base overlap, 1 base apart a zero-length DB
        let length = upper.position() - lower.position() - 1;
        let link_length = if lower.forward_orient() { length } else { -length };
        DbPair { lower, upper, link_length }
    }
    pub fn lower(&self) -> &Be {
        &self.lower
    }
    pub fn upper(&self) -> &Be {
        &self.upper
    }
    pub fn lower_sv(&self) -> &Sv {
        self.lower.sv()
    }
    pub fn upper_sv(&self) -> &Sv {
        self.upper.sv()
    }
    pub fn lower_link_on_start(&self) -> bool {
        self.lower.uses_start()
    }
    pub fn upper_link_on_start(&self) -> bool {
        self.upper.uses_start()
    }
    pub fn breakend(&self, is_start: bool) -> &Be {
        if is_start { &self.lower } else { &self.upper }
    }
    pub fn chromosome(&self) -> &str {
        self.lower.chromosome()
    }
    pub fn length(&self) -> i32 {
        self.link_length
    }
    pub fn has_breakend_of(&self, var: &Sv, use_start: bool) -> bool {
        (var == self.lower.sv() && self.lower.uses_start() == use_start) || (var == self.upper.sv() && self.upper.uses_start() == use_start)
    }
    pub fn has_breakend(&self, breakend: &Be) -> bool {
        breakend == &self.lower || breakend == &self.upper
    }
    pub fn other_sv(&self, var: &Sv) -> Sv {
        if self.lower.sv() == var { self.upper.sv().clone() } else { self.lower.sv().clone() }
    }
    pub fn other_breakend(&self, breakend: &Be) -> Be {
        if &self.lower == breakend { self.upper.clone() } else { self.lower.clone() }
    }
}

impl std::fmt::Display for DbPair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {}:{}:{} & {} {}:{}:{}",
            self.lower.sv().id(),
            self.lower.chromosome(),
            self.lower.position(),
            if self.lower.uses_start() { "start" } else { "end" },
            self.upper.sv().id(),
            self.upper.chromosome(),
            self.upper.position(),
            if self.upper.uses_start() { "start" } else { "end" }
        )
    }
}

/// ArmGroup: the SVs of a cluster on one chromosomal arm
pub struct ArmGroup {
    pub id: String,
    pub chromosome: String,
    pub arm: Arm,
    pub svs: Vec<Sv>,
    pub breakends: Vec<Be>,
    pub start_pos: i32,
    pub end_pos: i32,
}

impl ArmGroup {
    pub fn new(chr: &str, arm: Arm) -> ArmGroup {
        ArmGroup { id: make_chr_arm_str(chr, arm), chromosome: chr.to_string(), arm, svs: Vec::new(), breakends: Vec::new(), start_pos: -1, end_pos: -1 }
    }
    pub fn add_variant(&mut self, var: &Sv) {
        self.svs.push(var.clone());
        for be in [SE_START, SE_END] {
            if var.is_sgl_breakend() && be == SE_END {
                continue;
            }
            let breakend = var.breakend_se(be).unwrap();
            if breakend.chromosome() != self.chromosome || breakend.arm() != self.arm {
                continue;
            }
            let mut index = 0;
            while index < self.breakends.len() {
                if breakend.position() < self.breakends[index].position() {
                    break;
                }
                index += 1;
            }
            self.breakends.insert(index, breakend);
        }
        self.start_pos = self.breakends[0].position();
        self.end_pos = self.breakends[self.breakends.len() - 1].position();
    }
    pub fn matches(&self, other: &ArmGroup) -> bool {
        self.chromosome == other.chromosome && self.arm == other.arm
    }
}

/// ArmCluster: proximate breakends of a cluster on one arm
pub struct ArmCluster {
    pub id: i32,
    pub cluster_id: i32,
    pub breakends: Vec<Be>,
    pub chromosome: String,
    pub arm: Arm,
    pub start_pos: i32,
    pub end_pos: i32,
    pub arm_cluster_type: ArmClusterType,
    pub ti_count: i32,
}

impl ArmCluster {
    pub fn new(id: i32, cluster_id: i32, chr: &str, arm: Arm) -> ArmCluster {
        ArmCluster {
            id,
            cluster_id,
            breakends: Vec::new(),
            chromosome: chr.to_string(),
            arm,
            start_pos: 0,
            end_pos: 0,
            arm_cluster_type: ArmClusterType::UNKNOWN,
            ti_count: 0,
        }
    }
    pub fn chr_arm(&self) -> String {
        make_chr_arm_str(&self.chromosome, self.arm)
    }
    pub fn add_breakend(&mut self, breakend: &Be) {
        let mut index = 0;
        while index < self.breakends.len() {
            if breakend.position() < self.breakends[index].position() {
                break;
            }
            index += 1;
        }
        self.breakends.insert(index, breakend.clone());
        self.start_pos = self.breakends[0].position();
        self.end_pos = self.breakends[self.breakends.len() - 1].position();
    }
    pub fn type_str(&self) -> String {
        self.arm_cluster_type.to_string()
    }
}

impl std::fmt::Display for ArmCluster {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {} {}_{} {}:{}", self.id, self.arm_cluster_type, self.chromosome, self.arm, self.start_pos, self.end_pos)
    }
}

/// SvCluster
pub struct SvCluster {
    id: Cell<i32>,
    consistency_count: Cell<i32>,
    is_consistent: Cell<bool>,
    desc: RefCell<String>,
    type_counts: Cell<[i32; 7]>,
    is_resolved: Cell<bool>,
    resolved_type: Cell<ResolvedType>,
    annotation_list: RefCell<Vec<String>>,
    svs: RefCell<Vec<Sv>>,
    chains: RefCell<Vec<Chain>>,
    linked_pairs: RefCell<Vec<Pair>>,
    assembly_linked_pairs: RefCell<Vec<Pair>>,
    arm_groups: RefCell<Vec<ArmGroup>>,
    pub arm_clusters: RefCell<Vec<ArmCluster>>,
    chr_breakend_map: RefCell<JavaHashMap<String, Vec<Be>>>,
    unchained_svs: RefCell<Vec<Sv>>,
    loh_events: RefCell<Vec<crate::cn::Loh>>,
    clustering_reasons: RefCell<String>,
    requires_replication: Cell<bool>,
    long_del_dups: RefCell<Vec<Sv>>,
    foldbacks: RefCell<Vec<Sv>>,
    inversions: RefCell<Vec<Sv>>,
    double_minute_svs: RefCell<Vec<Sv>>,
    double_minute_chains: RefCell<Vec<Chain>>,
    has_linking_line_elements: Cell<bool>,
    requires_recalc: Cell<bool>,
    recalc_remote_sv_status: Cell<bool>,
    short_ti_remote_svs: RefCell<Vec<Sv>>,
    unlinked_remote_svs: RefCell<Vec<Sv>>,
    jcn_range: Cell<[f64; 2]>,
    pub metrics: RefCell<crate::analysis::ClusterMetrics>,
}

pub const CLUSTER_ANNOT_DM: &str = "DM";
pub const CLUSTER_ANNOT_BFB: &str = "BFB";
pub const CLUSTER_ANNOT_REP_REPAIR: &str = "REP_REPAIR";
pub const CLUSTER_ANNOT_SHATTERING: &str = "CHROMO_THRIP";

impl SvCluster {
    pub fn id(&self) -> i32 {
        self.id.get()
    }
    pub fn sv_count(&self) -> usize {
        self.svs.borrow().len()
    }
    pub fn desc(&self) -> String {
        self.desc.borrow().clone()
    }
    pub fn svs(&self) -> Vec<Sv> {
        self.svs.borrow().clone()
    }
    pub fn svs_ref(&self) -> Ref<'_, Vec<Sv>> {
        self.svs.borrow()
    }
    pub fn sv(&self, index: usize) -> Option<Sv> {
        self.svs.borrow().get(index).cloned()
    }
    pub fn sv0(&self) -> Sv {
        self.svs.borrow()[0].clone()
    }
    pub fn arm_groups(&self) -> Ref<'_, Vec<ArmGroup>> {
        self.arm_groups.borrow()
    }
    pub fn chr_breakend_map(&self) -> Ref<'_, JavaHashMap<String, Vec<Be>>> {
        self.chr_breakend_map.borrow()
    }
    /// a copy of the chromosome breakend map's entries in its iteration order
    pub fn chr_breakend_entries(&self) -> Vec<(String, Vec<Be>)> {
        self.chr_breakend_map.borrow().iter().map(|(k, v)| (k.clone(), v.clone())).collect()
    }
    pub fn requires_replication(&self) -> bool {
        self.requires_replication.get()
    }
    pub fn set_requires_replication(&self) {
        self.requires_replication.set(true)
    }
    pub fn add_loh_event(&self, e: &crate::cn::Loh) {
        let mut l = self.loh_events.borrow_mut();
        if !contains(&l, e) {
            l.push(e.clone());
        }
    }
    pub fn loh_events(&self) -> Vec<crate::cn::Loh> {
        self.loh_events.borrow().clone()
    }
    pub fn chains(&self) -> Vec<Chain> {
        self.chains.borrow().clone()
    }
    pub fn chain_count(&self) -> usize {
        self.chains.borrow().len()
    }
    pub fn chains_mut(&self) -> RefMut<'_, Vec<Chain>> {
        self.chains.borrow_mut()
    }

    pub fn add_chain(&self, chain: &Chain, reset_id: bool) {
        if reset_id {
            chain.set_id(self.chains.borrow().len() as i32);
        }
        self.chains.borrow_mut().push(chain.clone());
        let chain_svs = chain.sv_list();
        let mut u = self.unchained_svs.borrow_mut();
        for var in &chain_svs {
            remove_item(&mut u, var);
        }
    }

    pub fn is_fully_chained(&self, require_consistency: bool) -> bool {
        if !self.unchained_svs.borrow().is_empty() || self.chains.borrow().is_empty() {
            return false;
        }
        if require_consistency {
            for chain in self.chains.borrow().iter() {
                if !chain.is_consistent_check(require_consistency) {
                    return false;
                }
            }
        }
        true
    }

    pub fn dissolve_links_and_chains(&self) {
        let svs = self.svs.borrow().clone();
        *self.unchained_svs.borrow_mut() = svs;
        self.chains.borrow_mut().clear();
    }

    pub fn unlinked_svs(&self) -> Vec<Sv> {
        self.unchained_svs.borrow().clone()
    }
    pub fn unlinked_sv_count(&self) -> usize {
        self.unchained_svs.borrow().len()
    }
    pub fn linked_pairs(&self) -> Vec<Pair> {
        self.linked_pairs.borrow().clone()
    }
    pub fn assembly_linked_pairs(&self) -> Vec<Pair> {
        self.assembly_linked_pairs.borrow().clone()
    }
    pub fn set_assembly_linked_pairs(&self, pairs: Vec<Pair>) {
        *self.assembly_linked_pairs.borrow_mut() = pairs;
    }

    pub fn add_cluster_reason(&self, reason: ClusteringReason) {
        let name = reason.name();
        let mut r = self.clustering_reasons.borrow_mut();
        if !r.contains(&name) {
            crate::utils::append_str(&mut r, &name, ';');
        }
    }
    pub fn clustering_reasons(&self) -> String {
        self.clustering_reasons.borrow().clone()
    }
    pub fn has_cluster_reason(&self, reason: ClusteringReason) -> bool {
        self.clustering_reasons.borrow().contains(&reason.name())
    }

    pub fn is_consistent(&self) -> bool {
        self.update_cluster_details();
        self.is_consistent.get()
    }
    pub fn consistency_count(&self) -> i32 {
        self.update_cluster_details();
        self.consistency_count.get()
    }

    pub fn set_resolved(&self, toggle: bool, t: ResolvedType) {
        self.is_resolved.set(toggle);
        self.resolved_type.set(t);
        if self.desc.borrow().is_empty() {
            let d = self.cluster_types_as_string();
            *self.desc.borrow_mut() = d;
        }
    }
    pub fn is_resolved(&self) -> bool {
        self.is_resolved.get()
    }
    pub fn resolved_type(&self) -> ResolvedType {
        self.resolved_type.get()
    }

    pub fn update_cluster_details(&self) {
        if !self.requires_recalc.get() {
            return;
        }
        let svs = self.svs.borrow().clone();
        let c = crate::utils::calc_consistency_list(&svs);
        self.consistency_count.set(c);
        self.is_consistent.set(c == 0);
        let d = self.cluster_types_as_string();
        *self.desc.borrow_mut() = d;
        self.set_jcn_range();
        self.reset_breakend_map_indices();
        self.requires_recalc.set(false);
    }

    pub fn arm_count(&self) -> usize {
        self.arm_groups.borrow().len()
    }

    pub fn cluster_types_as_string(&self) -> String {
        let svs = self.svs.borrow();
        if svs.len() == 1 {
            return svs[0].type_str().to_string();
        }
        crate::utils::get_sv_types_str(&self.type_counts.get())
    }

    pub fn sgl_breakend_count(&self) -> i32 {
        let t = self.type_counts.get();
        t[SGL.ordinal()] + t[INF.ordinal()]
    }
    pub fn type_count(&self, t: StructuralVariantType) -> i32 {
        self.type_counts.get()[t.ordinal()]
    }

    pub fn long_del_dups(&self) -> Vec<Sv> {
        self.long_del_dups.borrow().clone()
    }
    pub fn foldbacks(&self) -> Vec<Sv> {
        self.foldbacks.borrow().clone()
    }
    pub fn foldback_count(&self) -> usize {
        self.foldbacks.borrow().len()
    }
    pub fn inversions(&self) -> Vec<Sv> {
        self.inversions.borrow().clone()
    }
    pub fn register_foldback(&self, var: &Sv) {
        let mut f = self.foldbacks.borrow_mut();
        if !contains(&f, var) {
            f.push(var.clone());
        }
    }
    pub fn deregister_foldback(&self, var: &Sv) {
        remove_item(&mut self.foldbacks.borrow_mut(), var);
    }
    pub fn register_inversion(&self, var: &Sv) {
        let mut f = self.inversions.borrow_mut();
        if !contains(&f, var) {
            f.push(var.clone());
        }
    }
    pub fn register_long_del_dup(&self, var: &Sv) {
        let mut f = self.long_del_dups.borrow_mut();
        if !contains(&f, var) {
            f.push(var.clone());
        }
    }

    pub fn double_minute_svs(&self) -> Vec<Sv> {
        self.double_minute_svs.borrow().clone()
    }
    pub fn double_minute_chains(&self) -> Vec<Chain> {
        self.double_minute_chains.borrow().clone()
    }
    pub fn set_double_minute_data(&self, svs: Vec<Sv>, chains: Vec<Chain>) {
        *self.double_minute_svs.borrow_mut() = svs;
        *self.double_minute_chains.borrow_mut() = chains;
    }
    pub fn clear_double_minute_data(&self) {
        self.double_minute_svs.borrow_mut().clear();
        self.double_minute_chains.borrow_mut().clear();
        remove_item(&mut self.annotation_list.borrow_mut(), &CLUSTER_ANNOT_DM.to_string());
    }

    pub fn mark_as_line(&self) {
        self.has_linking_line_elements.set(true);
        self.set_resolved(true, ResolvedType::LINE);
    }
    pub fn has_linking_line_elements(&self) -> bool {
        self.has_linking_line_elements.get()
    }

    pub fn unlinked_remote_svs(&self) -> Vec<Sv> {
        self.unlinked_remote_svs.borrow().clone()
    }

    pub fn set_arm_links(&self) {
        if !self.recalc_remote_sv_status.get() {
            return;
        }
        let mut short_ti: Vec<Sv> = Vec::new();
        for chain in self.chains.borrow().iter() {
            for pair in chain.linked_pairs().iter() {
                if pair.first().is_cross_arm() && pair.second().is_cross_arm() && pair.base_length() <= SHORT_TI_LENGTH {
                    if !contains(&short_ti, &pair.first()) {
                        short_ti.push(pair.first());
                    }
                    if !contains(&short_ti, &pair.second()) {
                        short_ti.push(pair.second());
                    }
                }
            }
        }
        let unlinked: Vec<Sv> =
            self.svs.borrow().iter().filter(|x| x.is_cross_arm()).filter(|x| !x.in_line_element()).filter(|x| !contains(&short_ti, x)).cloned().collect();
        *self.short_ti_remote_svs.borrow_mut() = short_ti;
        *self.unlinked_remote_svs.borrow_mut() = unlinked;
        self.recalc_remote_sv_status.set(false);
    }

    fn set_jcn_range(&self) {
        let svs = self.svs.borrow();
        let mut range = [0.0f64; 2];
        if svs.len() == 1 {
            let j = svs[0].jcn();
            self.jcn_range.set([j, j]);
            return;
        }
        range[SE_START] = -1.0;
        range[SE_END] = 0.0;
        let mut sv_calc_jcn_count = 0;
        let mut jcn_frequency: JavaHashMap<i32, i32> = JavaHashMap::new();
        let mut tightest_min_jcn = 0.0f64;
        let mut tightest_max_jcn = -1.0f64;
        let mut count_half_to_one_jcn = 0;
        let mut min_sv_jcn = -1.0f64;
        let mut max_assembled_multiple = 1;
        for var in svs.iter() {
            let sv_jcn = var.implied_jcn();
            max_assembled_multiple = max_assembled_multiple.max(var.max_assembled_breakend());
            if range[SE_START] < 0.0 || (sv_jcn as f64) < range[SE_START] {
                range[SE_START] = sv_jcn as f64;
                min_sv_jcn = sv_jcn as f64;
            }
            range[SE_END] = range[SE_END].max(sv_jcn as f64);
            if var.has_calculated_jcn() {
                sv_calc_jcn_count += 1;
                let min_jcn_int = var.jcn_min().ceil() as i32;
                let mut max_jcn_int = var.jcn_max().floor() as i32;
                max_jcn_int = min_jcn_int.max(max_jcn_int);
                if tightest_max_jcn == -1.0 || var.jcn_max() < tightest_max_jcn {
                    tightest_max_jcn = var.jcn_max();
                }
                tightest_min_jcn = var.jcn_min().max(tightest_min_jcn);
                if var.jcn_min() < 1.0 && var.jcn_max() > 0.5 {
                    count_half_to_one_jcn += 1;
                }
                for i in min_jcn_int..=max_jcn_int {
                    let c = jcn_frequency.get(&i).copied();
                    jcn_frequency.insert(i, c.map(|c| c + 1).unwrap_or(1));
                }
            }
        }
        if sv_calc_jcn_count > 0 {
            range[SE_START] = -1.0;
            range[SE_END] = 0.0;
            for (ploidy, sv_count) in jcn_frequency.iter() {
                let ploidy = *ploidy as f64;
                if *sv_count == sv_calc_jcn_count {
                    range[SE_END] = ploidy;
                    range[SE_START] = ploidy;
                    break;
                }
                if ploidy > 0.0 && (range[SE_START] < 0.0 || ploidy < range[SE_START]) {
                    range[SE_START] = ploidy;
                }
                range[SE_START] = range[SE_START].max(min_sv_jcn);
                range[SE_END] = range[SE_END].max(ploidy);
            }
            if range[SE_START] < range[SE_END] && max_assembled_multiple == 1 {
                if tightest_max_jcn > tightest_min_jcn && tightest_max_jcn - tightest_min_jcn < 1.0 {
                    range = [1.0, 1.0];
                } else if count_half_to_one_jcn == sv_calc_jcn_count {
                    range = [1.0, 1.0];
                }
            }
        }
        if max_assembled_multiple > 1 {
            range[SE_END] = range[SE_END].max(max_assembled_multiple as f64 * range[SE_START]);
        }
        self.jcn_range.set(range);
    }

    pub fn jcn_range(&self) -> [f64; 2] {
        self.jcn_range.get()
    }
    pub fn min_jcn(&self) -> f64 {
        let svs = self.svs.borrow();
        if svs.is_empty() {
            return 0.0;
        }
        svs.iter().map(|x| x.jcn()).fold(f64::INFINITY, crate::utils::jmin)
    }
    pub fn max_jcn(&self) -> f64 {
        let svs = self.svs.borrow();
        if svs.is_empty() {
            return 0.0;
        }
        svs.iter().map(|x| x.jcn()).fold(f64::NEG_INFINITY, crate::utils::jmax)
    }

    pub fn has_varied_jcn(&self) -> bool {
        if self.requires_recalc.get() {
            self.update_cluster_details();
        }
        if self.svs.borrow().len() == 1 {
            return false;
        }
        let r = self.jcn_range.get();
        r[SE_END] > r[SE_START] && r[SE_START] >= 0.0
    }

    fn reset_breakend_map_indices(&self) {
        for (_, list) in self.chr_breakend_map.borrow().iter() {
            for (i, b) in list.iter().enumerate() {
                b.set_cluster_chr_pos_index(i as i32);
            }
        }
    }

    pub fn cache_linked_pairs(&self) {
        let mut pairs: Vec<Pair> = Vec::new();
        for chain in self.chains.borrow().iter() {
            for pair in chain.linked_pairs() {
                if pairs.iter().any(|x| x.matches(&pair)) {
                    continue;
                }
                pairs.push(pair.clone());
                pair.first().add_linked_pair(&pair, pair.first_link_on_start());
                pair.second().add_linked_pair(&pair, pair.second_link_on_start());
            }
        }
        *self.linked_pairs.borrow_mut() = pairs;
    }

    pub fn find_chain(&self, var: &Sv) -> Option<Chain> {
        self.chains.borrow().iter().find(|c| c.sv_index(var) >= 0).cloned()
    }
    pub fn find_chains(&self, var: &Sv) -> Vec<Chain> {
        self.chains.borrow().iter().filter(|c| c.has_sv(var)).cloned().collect()
    }
    pub fn find_same_chain_for_svs(&self, var1: &Sv, var2: &Sv) -> Option<Chain> {
        let c1 = self.find_chains(var1);
        let c2 = self.find_chains(var2);
        for a in &c1 {
            for b in &c2 {
                if a == b {
                    return Some(a.clone());
                }
            }
        }
        None
    }

    pub fn chain_id(&self, var: &Sv) -> i32 {
        if let Some(chain) = self.find_chain(var) {
            return chain.id();
        }
        let u = self.unchained_svs.borrow();
        for (i, uv) in u.iter().enumerate() {
            if uv == var {
                return (self.chains.borrow().len() + i + 1) as i32;
            }
        }
        var.id()
    }

    pub fn annotation_list(&self) -> Vec<String> {
        self.annotation_list.borrow().clone()
    }
    pub fn add_annotation(&self, a: &str) {
        let mut l = self.annotation_list.borrow_mut();
        if l.iter().any(|x| x == a) {
            return;
        }
        l.push(a.to_string());
    }
    pub fn has_annotation(&self, a: &str) -> bool {
        self.annotation_list.borrow().iter().any(|x| x == a)
    }
    pub fn annotations(&self) -> String {
        self.annotation_list.borrow().join(";")
    }

    pub fn determine_requires_replication(&self) {
        if !self.has_varied_jcn() && !self.requires_replication() {
            return;
        }
        let cluster_min_jcn = self.jcn_range.get()[SE_START];
        if cluster_min_jcn <= 0.0 {
            hmf_common::debug!(
                "cluster({}) warning: invalid JCN variation(min={} max={})",
                self.id(),
                hmf_common::fmt::java_double_to_string(cluster_min_jcn),
                hmf_common::fmt::java_double_to_string(self.jcn_range.get()[SE_END])
            );
            return;
        }
        for var in self.svs.borrow().iter() {
            let sv_jcn = var.implied_jcn();
            let max_assembly_breakends = var.max_assembled_breakend();
            let mut sv_multiple = crate::utils::jround(sv_jcn as f64 / cluster_min_jcn) as i32;
            if max_assembly_breakends > 1 {
                sv_multiple = sv_multiple.max(max_assembly_breakends);
            }
            if sv_multiple <= 1 {
                continue;
            }
            if !self.requires_replication() {
                self.set_requires_replication();
            }
        }
    }

    pub fn is_synthetic_type(&self) -> bool {
        if self.svs.borrow().len() == 1 {
            return false;
        }
        crate::analysis::classification::is_synthetic_type(self)
    }
}



impl Cluster {
    pub fn create(cluster_id: i32) -> Cluster {
        Cluster::new(SvCluster {
            id: Cell::new(cluster_id),
            consistency_count: Cell::new(0),
            is_consistent: Cell::new(false),
            desc: RefCell::new(String::new()),
            type_counts: Cell::new([0; 7]),
            is_resolved: Cell::new(false),
            resolved_type: Cell::new(ResolvedType::NONE),
            annotation_list: RefCell::new(Vec::new()),
            svs: RefCell::new(Vec::new()),
            chains: RefCell::new(Vec::new()),
            linked_pairs: RefCell::new(Vec::new()),
            assembly_linked_pairs: RefCell::new(Vec::new()),
            arm_groups: RefCell::new(Vec::new()),
            arm_clusters: RefCell::new(Vec::new()),
            chr_breakend_map: RefCell::new(JavaHashMap::new()),
            unchained_svs: RefCell::new(Vec::new()),
            loh_events: RefCell::new(Vec::new()),
            clustering_reasons: RefCell::new(String::new()),
            requires_replication: Cell::new(false),
            long_del_dups: RefCell::new(Vec::new()),
            foldbacks: RefCell::new(Vec::new()),
            inversions: RefCell::new(Vec::new()),
            double_minute_svs: RefCell::new(Vec::new()),
            double_minute_chains: RefCell::new(Vec::new()),
            has_linking_line_elements: Cell::new(false),
            requires_recalc: Cell::new(true),
            recalc_remote_sv_status: Cell::new(false),
            short_ti_remote_svs: RefCell::new(Vec::new()),
            unlinked_remote_svs: RefCell::new(Vec::new()),
            jcn_range: Cell::new([0.0, 0.0]),
            metrics: RefCell::new(crate::analysis::ClusterMetrics::default()),
        })
    }

    pub fn add_variant(&self, var: &Sv) {
        if contains(&self.svs.borrow(), var) {
            hmf_common::error!("cluster({}) attempting to add SV({}) again", self.id(), var.id());
            return;
        }
        self.svs.borrow_mut().push(var.clone());
        self.requires_recalc.set(true);
        self.annotation_list.borrow_mut().clear();
        self.double_minute_svs.borrow_mut().clear();
        self.double_minute_chains.borrow_mut().clear();
        var.set_cluster(Some(self.clone()));
        self.unchained_svs.borrow_mut().push(var.clone());
        if !self.has_linking_line_elements.get() {
            self.is_resolved.set(false);
            self.resolved_type.set(ResolvedType::NONE);
        }
        let mut tc = self.type_counts.get();
        tc[var.sv_type().ordinal()] += 1;
        self.type_counts.set(tc);
        if var.sv_type() == BND || var.is_cross_arm() {
            self.recalc_remote_sv_status.set(true);
        }
        crate::utils::add_sv_to_chr_breakend_map(var, &mut self.chr_breakend_map.borrow_mut());
        for be in [SE_START, SE_END] {
            if be == SE_END && var.is_sgl_breakend() {
                continue;
            }
            if be == SE_END && var.is_local() {
                continue;
            }
            let use_start = is_start(be);
            let mut groups = self.arm_groups.borrow_mut();
            let mut found = false;
            for g in groups.iter_mut() {
                if g.chromosome == var.chromosome(use_start) && g.arm == var.arm(use_start) {
                    g.add_variant(var);
                    found = true;
                    break;
                }
            }
            if !found {
                let mut g = ArmGroup::new(var.chromosome(use_start), var.arm(use_start));
                g.add_variant(var);
                groups.push(g);
            }
        }
    }

    pub fn merge_other_cluster(&self, other: &Cluster) {
        self.merge_other_cluster_log(other, true);
    }

    pub fn merge_other_cluster_log(&self, other: &Cluster, log_details: bool) {
        if other == self || other.id() == self.id() {
            hmf_common::error!("attempting to merge same cluster({})", self.id());
            return;
        }
        if other.sv_count() > self.sv_count() || other.has_annotation(CLUSTER_ANNOT_DM) {
            if log_details {
                hmf_common::debug!(
                    "cluster({} svs={}) merges in other cluster({} svs={}) and adopts new ID",
                    self.id(),
                    self.sv_count(),
                    other.id(),
                    other.sv_count()
                );
            }
            self.id.set(other.id());
        } else if log_details {
            hmf_common::debug!("cluster({} svs={}) merges in other cluster({} svs={})", self.id(), self.sv_count(), other.id(), other.sv_count());
        }
        self.add_variant_lists(other);
        for loh in other.loh_events() {
            self.add_loh_event(&loh);
        }
        let reasons = other.clustering_reasons();
        for r in reasons.split(';') {
            if r.is_empty() {
                break;
            }
            self.add_cluster_reason(ClusteringReason::value_of(r).expect("clustering reason"));
        }
        if other.has_linking_line_elements() {
            self.mark_as_line();
        }
        let other_chains = other.chains();
        for c in other_chains.iter().filter(|c| c.assembly_link_count() == c.link_count()) {
            self.chains.borrow_mut().push(c.clone());
        }
    }

    fn add_variant_lists(&self, other: &Cluster) {
        for var in other.svs() {
            self.add_variant(&var);
        }
        let ap = other.assembly_linked_pairs();
        self.assembly_linked_pairs.borrow_mut().extend(ap);
        if other.requires_replication() {
            self.requires_replication.set(true);
        }
        let inv = other.inversions();
        self.inversions.borrow_mut().extend(inv);
        let fb = other.foldbacks();
        self.foldbacks.borrow_mut().extend(fb);
        let ldd = other.long_del_dups();
        self.long_del_dups.borrow_mut().extend(ldd);
    }

    pub fn find_arm_cluster_id(&self, breakend: &Be) -> Option<usize> {
        self.arm_clusters.borrow().iter().position(|ac| contains(&ac.breakends, breakend))
    }
}

impl std::fmt::Display for SvCluster {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: SVs({}) res({}) desc({}) chains({})", self.id(), self.sv_count(), self.resolved_type(), self.desc(), self.chain_count())
    }
}

/// SvChain: linked pairs in order, each pair's 'first' SV linking to the preceding pair
pub struct SvChain {
    id: Cell<i32>,
    sv_list: RefCell<Vec<Sv>>,
    linked_pairs: RefCell<Vec<Pair>>,
    open_breakends: RefCell<[Option<Be>; 2]>,
    jcn: Cell<f64>,
    jcn_uncertainty: Cell<f64>,
    is_closed_loop: Cell<bool>,
    double_minute: Cell<bool>,
    link_sum: Cell<i32>,
}

impl SvChain {
    pub fn id(&self) -> i32 {
        self.id.get()
    }
    pub fn set_id(&self, id: i32) {
        self.id.set(id)
    }
    pub fn sv_count(&self) -> usize {
        self.sv_list.borrow().len()
    }
    pub fn sv_list(&self) -> Vec<Sv> {
        self.sv_list.borrow().clone()
    }
    pub fn link_count(&self) -> usize {
        self.linked_pairs.borrow().len()
    }
    pub fn linked_pairs(&self) -> Vec<Pair> {
        self.linked_pairs.borrow().clone()
    }
    pub fn linked_pairs_ref(&self) -> Ref<'_, Vec<Pair>> {
        self.linked_pairs.borrow()
    }
    pub fn linked_pair(&self, i: usize) -> Pair {
        self.linked_pairs.borrow()[i].clone()
    }
    /// getLinkedPairs() for callers that edit the list directly
    pub fn linked_pairs_mut(&self) -> RefMut<'_, Vec<Pair>> {
        self.linked_pairs.borrow_mut()
    }

    pub fn add_link(&self, pair: &Pair, add_to_start: bool) {
        {
            let mut l = self.sv_list.borrow_mut();
            if !contains(&l, &pair.first()) {
                l.push(pair.first());
            }
            if !contains(&l, &pair.second()) {
                l.push(pair.second());
            }
        }
        let new_pair = Pair::copy_of(pair);
        self.link_sum.set(self.link_sum.get() + pair.first().id() + pair.second().id());
        if self.linked_pairs.borrow().is_empty() {
            self.linked_pairs.borrow_mut().push(new_pair);
            self.update_breakends();
            return;
        }
        if (add_to_start && !self.can_add_linked_pair_to_start(&new_pair)) || (!add_to_start && !self.can_add_linked_pair_to_end(&new_pair)) {
            hmf_common::error!("chain({}) trying to add invalid link({}) to {}", self.id(), new_pair, if add_to_start { "start" } else { "end" });
            return;
        }
        let switch = {
            let lp = self.linked_pairs.borrow();
            (add_to_start && new_pair.second() != lp[0].first()) || (!add_to_start && new_pair.first() != lp[lp.len() - 1].second())
        };
        if switch {
            new_pair.switch_svs();
        }
        if add_to_start {
            self.linked_pairs.borrow_mut().insert(0, new_pair);
        } else {
            self.linked_pairs.borrow_mut().push(new_pair);
        }
        self.update_breakends();
    }

    pub fn set_jcn_data(&self, jcn: f64, uncertainty: f64) {
        self.jcn.set(jcn);
        self.jcn_uncertainty.set(uncertainty);
    }
    pub fn jcn(&self) -> f64 {
        self.jcn.get()
    }
    pub fn jcn_uncertainty(&self) -> f64 {
        self.jcn_uncertainty.get()
    }

    pub fn chain_end_sv(&self, is_start: bool) -> Option<Sv> {
        let lp = self.linked_pairs.borrow();
        if lp.is_empty() {
            return None;
        }
        if is_start { Some(lp[0].first()) } else { Some(lp[lp.len() - 1].second()) }
    }
    pub fn first_sv(&self) -> Sv {
        self.chain_end_sv(true).expect("chain first SV")
    }
    pub fn last_sv(&self) -> Sv {
        self.chain_end_sv(false).expect("chain last SV")
    }
    pub fn first_link_open_on_start(&self) -> bool {
        self.linked_pairs.borrow().first().map(|p| p.first_unlinked_on_start()).unwrap_or(false)
    }
    pub fn last_link_open_on_start(&self) -> bool {
        self.linked_pairs.borrow().last().map(|p| !p.second_link_on_start()).unwrap_or(false)
    }

    fn update_breakends(&self) {
        let s = self.first_sv().breakend(self.first_link_open_on_start());
        let e = self.last_sv().breakend(self.last_link_open_on_start());
        *self.open_breakends.borrow_mut() = [s, e];
    }

    pub fn open_breakend(&self, is_start: bool) -> Option<Be> {
        self.open_breakends.borrow()[se(is_start)].clone()
    }

    pub fn can_add_linked_pair_to_start(&self, pair: &Pair) -> bool {
        self.can_add_linked_pair(pair, true, false)
    }
    pub fn can_add_linked_pair_to_end(&self, pair: &Pair) -> bool {
        self.can_add_linked_pair(pair, false, false)
    }
    pub fn can_add_linked_pair(&self, pair: &Pair, to_start: bool, _allow_replicated: bool) -> bool {
        if self.linked_pairs.borrow().is_empty() {
            return true;
        }
        if contains(&self.linked_pairs.borrow(), pair) {
            return false;
        }
        let chain_sv = if to_start { self.first_sv() } else { self.last_sv() };
        let chain_open_side = if to_start { self.first_link_open_on_start() } else { self.last_link_open_on_start() };
        if pair.first() == chain_sv && pair.first_unlinked_on_start() != chain_open_side {
            true
        } else {
            pair.second() == chain_sv && pair.second_unlinked_on_start() != chain_open_side
        }
    }

    pub fn link_would_close_chain(&self, pair: &Pair) -> bool {
        let chain_start = self.open_breakend(true);
        let chain_end = self.open_breakend(false);
        let (f, s) = (Some(pair.first_breakend()), Some(pair.second_breakend()));
        (f == chain_start && s == chain_end) || (f == chain_end && s == chain_start)
    }

    pub fn could_close_chain(&self) -> bool {
        {
            let l = self.sv_list.borrow();
            if l.len() == 1 && l[0].sv_type() == DUP {
                return true;
            }
        }
        let chain_start = self.open_breakend(true);
        let chain_end = self.open_breakend(false);
        if let (Some(cs), Some(ce)) = (&chain_start, &chain_end) {
            if !cs.sv().is_sgl_breakend() && !ce.sv().is_sgl_breakend() {
                if ce.sv() == cs.sv() && cs != ce {
                    return true;
                }
                let min_ti_length = crate::chaining::utils::get_min_templated_insertion_length(cs, ce);
                return crate::chaining::utils::is_possible_link(
                    cs.chromosome(),
                    cs.position(),
                    cs.orientation(),
                    ce.chromosome(),
                    ce.position(),
                    ce.orientation(),
                    min_ti_length,
                );
            }
        }
        false
    }

    pub fn is_closed_loop(&self) -> bool {
        let l = self.sv_list.borrow();
        self.is_closed_loop.get() || (l.len() == 1 && l[0].sv_type() == DUP)
    }
    pub fn is_double_minute(&self) -> bool {
        self.double_minute.get()
    }
    pub fn set_double_minute(&self, t: bool) {
        self.double_minute.set(t)
    }

    pub fn close_chain_with(&self, reason: &str, link_index: i32) {
        let chain_start = self.open_breakend(true).expect("chain start");
        let chain_end = self.open_breakend(false).expect("chain end");
        if chain_start.chromosome() != chain_end.chromosome() {
            hmf_common::error!("chain({}) cannot close chain with ends: {} & {}", self.id(), chain_start, chain_end);
            return;
        }
        let link = Pair::from_breakends(&chain_end, &chain_start);
        link.set_link_reason(reason, link_index);
        self.linked_pairs.borrow_mut().push(link);
        self.is_closed_loop.set(true);
        self.link_sum.set(self.link_sum.get() + chain_end.sv().id() + chain_start.sv().id());
    }
    pub fn close_chain(&self) {
        self.is_closed_loop.set(true);
    }

    pub fn is_consistent(&self) -> bool {
        self.is_consistent_check(false)
    }
    pub fn is_consistent_check(&self, check_centromere: bool) -> bool {
        if self.first_sv().is_sgl_breakend() || self.last_sv().is_sgl_breakend() {
            return false;
        }
        let ob = self.open_breakends.borrow();
        let (s, e) = (ob[0].as_ref().unwrap(), ob[1].as_ref().unwrap());
        if check_centromere {
            let mut traverses_centromere = false;
            for pair in self.linked_pairs.borrow().iter() {
                if pair.first_breakend().arm() != pair.second_breakend().arm() {
                    if traverses_centromere {
                        return false;
                    }
                    traverses_centromere = true;
                }
            }
            if traverses_centromere {
                return true;
            }
            if s.chr_arm() != e.chr_arm() {
                return false;
            }
            s.orientation() != e.orientation()
        } else {
            crate::utils::calc_consistency_be(s) + crate::utils::calc_consistency_be(e) == 0
        }
    }

    pub fn copy_from(&self, other: &SvChain) {
        for pair in other.linked_pairs() {
            self.add_link(&pair, false);
        }
    }

    pub fn has_sv(&self, var: &Sv) -> bool {
        contains(&self.sv_list.borrow(), var)
    }
    pub fn sv_index(&self, var: &Sv) -> i32 {
        index_of(&self.sv_list.borrow(), var)
    }

    pub fn sv_indices(&self, var: &Sv) -> String {
        let mut s = String::new();
        for (index, pair) in self.linked_pairs.borrow().iter().enumerate() {
            let info = if &pair.first() == var {
                format!("{}{}", index, if pair.first_link_on_start() { "s" } else { "e" })
            } else if &pair.second() == var {
                format!("{}{}", index, if pair.second_link_on_start() { "s" } else { "e" })
            } else {
                continue;
            };
            crate::utils::append_str(&mut s, &info, ';');
        }
        s
    }

    pub fn length(&self, close_ends: bool) -> i64 {
        let lp = self.linked_pairs.borrow();
        let mut length: i64 = lp.iter().map(|x| (x.base_length() as i64).abs()).sum();
        if close_ends {
            if let (Some(cs), Some(ce)) = (self.open_breakend(true), self.open_breakend(false)) {
                if lp.len() == 1 && ce.sv() == cs.sv() {
                    return length;
                }
                if cs.chromosome() == ce.chromosome() {
                    length += (cs.position() - ce.position()).abs() as i64;
                }
            }
        }
        length
    }

    pub fn has_repeated_sv(&self) -> bool {
        let n = self.linked_pairs.borrow().len();
        let m = self.sv_list.borrow().len();
        if self.is_closed_loop.get() { n > m } else { n >= m }
    }
    pub fn assembly_link_count(&self) -> usize {
        self.linked_pairs.borrow().iter().filter(|x| x.is_assembled()).count()
    }
    pub fn has_matching_svs(&self, other: &SvChain) -> bool {
        let l = self.sv_list.borrow();
        other.sv_list().iter().all(|v| contains(&l, v))
    }
    pub fn same_links(&self, other: &SvChain) -> bool {
        if other.link_count() != self.link_count() {
            return false;
        }
        let o = other.linked_pairs();
        self.linked_pairs.borrow().iter().all(|p| o.iter().any(|x| x.matches(p)))
    }
    pub fn link_sum(&self) -> i32 {
        self.link_sum.get()
    }
}

impl Chain {
    pub fn create(chain_id: i32) -> Chain {
        Chain::new(SvChain {
            id: Cell::new(chain_id),
            sv_list: RefCell::new(Vec::new()),
            linked_pairs: RefCell::new(Vec::new()),
            open_breakends: RefCell::new([None, None]),
            jcn: Cell::new(0.0),
            jcn_uncertainty: Cell::new(0.0),
            is_closed_loop: Cell::new(false),
            double_minute: Cell::new(false),
            link_sum: Cell::new(0),
        })
    }
}

impl std::fmt::Display for SvChain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: svs({}) links({}) jcn({}) {}",
            self.id(),
            self.sv_count(),
            self.link_count(),
            hmf_common::fmt::java_fixed(self.jcn(), 1),
            if self.is_closed_loop.get() { "closed" } else { "" }
        )
    }
}

/// SvChain.checkIsValid
pub fn chain_is_valid(chain: &SvChain) -> bool {
    if chain.sv_count() == 0 || chain.link_count() == 0 {
        return false;
    }
    let lp = chain.linked_pairs();
    let svs = chain.sv_list();
    for i in 0..lp.len() - 1 {
        let pair = &lp[i];
        if pair.first() == pair.second() && !pair.is_dup_link() {
            return false;
        }
        if !contains(&svs, &pair.first()) || !contains(&svs, &pair.second()) {
            return false;
        }
        let next = &lp[i + 1];
        if pair.second() != next.first() {
            return false;
        }
        if pair.second_link_on_start() == next.first_link_on_start() {
            return false;
        }
    }
    true
}
