//! ChainJcnLimits and SegmentJcn: the allele and cluster JCN of each segment between a cluster's
//! breakends, used to limit the links a chain can make; and the JCN matching tests.

use hmf_common::java_hash::JavaHashMap;

use crate::cn::JcnCalcData;
use crate::types::*;
use crate::utils::{copy_numbers_equal, jmax, jmin, jround};

pub const EXHAUSTED_JCN_PERC: f64 = 0.1;
pub const EXHAUSTED_JCN_ABS: f64 = 0.2;

/// SegmentJcn
#[derive(Clone, Debug)]
pub struct SegmentJcn {
    pub pos_start: i32,
    pub pos_end: i32,
    pub major_ap: f64,
    pub minor_ap: f64,
    pub a_fixed_ap: f64,
    pub b_undisrupted_ap: f64,
    is_valid: bool,
    cluster_ap: f64,
    link_jcn: f64,
    exhausted_level: f64,
}

impl SegmentJcn {
    pub fn new(pos_start: i32, pos_end: i32, major: f64, minor: f64) -> Self {
        SegmentJcn {
            pos_start,
            pos_end,
            major_ap: major,
            minor_ap: minor,
            a_fixed_ap: 0.0,
            b_undisrupted_ap: 0.0,
            is_valid: false,
            cluster_ap: 0.0,
            link_jcn: 0.0,
            exhausted_level: 0.0,
        }
    }
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }
    pub fn set_valid(&mut self, t: bool) {
        self.is_valid = t
    }
    pub fn length(&self) -> i32 {
        self.pos_end - self.pos_start
    }
    pub fn cluster_jcn(&self) -> f64 {
        self.cluster_ap
    }
    pub fn set_cluster_jcn(&mut self, jcn: f64) {
        self.cluster_ap = jcn;
        self.exhausted_level = jmin(EXHAUSTED_JCN_PERC * jcn, EXHAUSTED_JCN_ABS);
    }
    pub fn unlinked_jcn(&self) -> f64 {
        jmax(self.cluster_ap - self.link_jcn, 0.0)
    }
    pub fn add_link_jcn(&mut self, link: f64) {
        if self.unlinked_jcn() > link {
            self.link_jcn += link;
        } else {
            self.link_jcn = self.cluster_ap;
        }
    }
}

pub const CLUSTER_ALLELE_JCN_MIN: f64 = 0.15;
pub const RANGE_TOTAL: usize = 0;
pub const DELETED_TOTAL: usize = 1;
const JCN_INTEGER_PROXIMITY: f64 = 0.25;

/// ChainJcnLimits
#[derive(Default)]
pub struct ChainJcnLimits {
    pub chr_allele_jcns: JavaHashMap<String, Vec<SegmentJcn>>,
    valid_allele_jcn_segment_perc: f64,
    chr_breakend_map: Vec<(String, Vec<Be>)>,
    cluster_id: i32,
}

impl ChainJcnLimits {
    pub fn new() -> Self {
        ChainJcnLimits { cluster_id: -1, ..Default::default() }
    }
    pub fn valid_allele_jcn_segment_perc(&self) -> f64 {
        self.valid_allele_jcn_segment_perc
    }
    /// the cluster's chromosome breakend map, as its entries in iteration order
    pub fn initialise(&mut self, cluster_id: i32, chr_breakend_map: Vec<(String, Vec<Be>)>) {
        self.chr_breakend_map = chr_breakend_map;
        self.cluster_id = cluster_id;
        self.chr_allele_jcns.clear();
        self.valid_allele_jcn_segment_perc = 0.0;
    }

    pub fn calc_range_data(&self) -> Option<[i64; 2]> {
        if self.chr_allele_jcns.is_empty() {
            return None;
        }
        let mut r = [0i64, 0i64];
        for (_, list) in self.chr_allele_jcns.iter() {
            for s in list {
                if s.cluster_jcn() >= CLUSTER_ALLELE_JCN_MIN {
                    r[RANGE_TOTAL] += s.length() as i64;
                } else {
                    r[DELETED_TOTAL] += s.length() as i64;
                }
            }
        }
        Some(r)
    }

    pub fn determine_breakend_jcns(&mut self) {
        let mut total_seg = 0;
        let mut total_valid = 0;
        let entries = std::mem::take(&mut self.chr_breakend_map);
        for (chromosome, list) in &entries {
            let mut allele: Vec<SegmentJcn> = Vec::with_capacity(list.len());
            let mut in_segment = false;
            let mut seg_start_breakend: Option<Be> = None;
            let mut seg_start_index = 0;
            for i in 0..list.len() {
                let breakend = &list[i];
                let next = list.get(i + 1);
                allele.push(SegmentJcn::new(
                    breakend.position(),
                    next.map(|n| n.position()).unwrap_or(breakend.position()),
                    breakend.major_allele_jcn(false),
                    breakend.minor_allele_jcn(false),
                ));
                if !in_segment {
                    in_segment = true;
                    seg_start_breakend = Some(breakend.clone());
                    seg_start_index = i;
                } else if next.is_none_or(|n| n.chr_pos_index() > breakend.chr_pos_index() + 1) {
                    in_segment = false;
                    let valid =
                        calculate_none_cluster_segment_ploidies(&mut allele, seg_start_index, i, seg_start_breakend.as_ref().unwrap(), breakend);
                    total_seg += 1;
                    if valid {
                        total_valid += 1;
                    }
                }
            }
            self.calculate_cluster_segment_jcns(&mut allele);
            self.chr_allele_jcns.insert(chromosome.clone(), allele);
        }
        hmf_common::debug!(
            "cluster({}) chromosomes({}) AP totalSegments({} valid={})",
            self.cluster_id,
            entries.len(),
            total_seg,
            total_valid
        );
        self.chr_breakend_map = entries;
    }

    fn calculate_cluster_segment_jcns(&mut self, allele: &mut [SegmentJcn]) {
        let b_min = allele.iter().filter(|x| x.is_valid()).map(|x| x.b_undisrupted_ap).fold(None, |acc: Option<f64>, v| {
            Some(match acc {
                None => v,
                Some(a) => jmin(a, v),
            })
        });
        let b_min = b_min.unwrap_or(0.0);
        let mut valid_segments = 0;
        let n = allele.len();
        for s in allele.iter_mut() {
            if !s.is_valid() {
                continue;
            }
            valid_segments += 1;
            let (major, minor, a) = (s.major_ap, s.minor_ap, s.a_fixed_ap);
            s.b_undisrupted_ap = b_min;
            let cluster_jcn = if major > a + 0.5 { major - b_min } else { minor - b_min };
            s.set_cluster_jcn(jmax(cluster_jcn, 0.0));
            s.set_valid(true);
        }
        self.valid_allele_jcn_segment_perc = valid_segments as f64 / n as f64;
    }

    /// linkHasJcnSupport
    pub fn link_has_jcn_support(&self, pair: &Pair, jcn: f64) -> bool {
        let lower = pair.breakend(true);
        let upper = pair.breakend(false);
        let Some(segments) = self.chr_allele_jcns.get(&pair.chromosome()) else { return false };
        if segments.is_empty() {
            return false;
        }
        let start = lower.cluster_chr_pos_index();
        for i in (start + 1)..upper.cluster_chr_pos_index() {
            let s = &segments[i as usize];
            if !s.is_valid() {
                continue;
            }
            let sj = s.unlinked_jcn();
            if sj < jcn && !copy_numbers_equal(sj, jcn) {
                return false;
            }
        }
        true
    }

    /// assignLinkJcn
    pub fn assign_link_jcn(&mut self, pair: &Pair, jcn: f64) {
        let lower = pair.breakend(true);
        let upper = pair.breakend(false);
        let Some(segments) = self.chr_allele_jcns.get_mut(&pair.chromosome()) else { return };
        if segments.is_empty() {
            return;
        }
        for i in lower.cluster_chr_pos_index()..upper.cluster_chr_pos_index() {
            let s = &mut segments[i as usize];
            if !s.is_valid() {
                continue;
            }
            s.add_link_jcn(jcn);
        }
    }
}

fn is_jcn_close_to_integer(jcn: f64) -> bool {
    let r = (jcn % 1.0).abs();
    r <= JCN_INTEGER_PROXIMITY || r >= (1.0 - JCN_INTEGER_PROXIMITY)
}

fn calculate_none_cluster_segment_ploidies(allele: &mut [SegmentJcn], start_index: usize, end_index: usize, start_breakend: &Be, end_breakend: &Be) -> bool {
    let start_major = start_breakend.major_allele_jcn(true);
    let start_minor = start_breakend.minor_allele_jcn(true);
    let start_major_r = jround(start_major) as i32;
    let start_minor_r = jround(start_minor) as i32;
    let start_major_int = is_jcn_close_to_integer(start_major);
    let start_minor_int = is_jcn_close_to_integer(start_minor);
    let mut freq: JavaHashMap<i32, i32> = JavaHashMap::new();
    let seg_count = end_index - start_index + 1;
    for s in &allele[start_index..=end_index] {
        let major = jround(s.major_ap) as i32;
        let minor = jround(s.minor_ap) as i32;
        let c = freq.get(&major).copied();
        freq.insert(major, c.map(|v| v + 1).unwrap_or(1));
        if minor != major {
            let c = freq.get(&minor).copied();
            freq.insert(minor, c.map(|v| v + 1).unwrap_or(1));
        }
    }
    let mut max_count = 0;
    let mut a_jcn = 0.0;
    for (k, v) in freq.iter() {
        if *v > max_count {
            max_count = *v;
            a_jcn = *k as f64;
        }
    }
    if (max_count as f64) < seg_count as f64 * 0.9 {
        return false;
    }
    let mut b_min = -1.0;
    let start_cluster_jcn = if start_breakend.forward_orient() { start_breakend.jcn() } else { 0.0 };
    if start_major_int && start_major_r as f64 == a_jcn {
        b_min = start_minor - start_cluster_jcn;
    } else if start_minor_int && start_minor_r as f64 == a_jcn {
        b_min = start_major - start_cluster_jcn;
    }
    let end_cluster_jcn = if end_breakend.reverse_orient() { end_breakend.jcn() } else { 0.0 };
    for i in start_index..=end_index {
        let mut major = jround(allele[i].major_ap) as i32;
        let mut minor = jround(allele[i].minor_ap) as i32;
        if major as f64 == a_jcn {
            if i == end_index {
                // int -= double: (int)(minor - endClusterJcn)
                minor = (minor as f64 - end_cluster_jcn) as i32;
            }
            b_min = if b_min == -1.0 { minor as f64 } else { jmin(minor as f64, b_min) };
        } else if minor as f64 == a_jcn {
            if i == end_index {
                major = (major as f64 - end_cluster_jcn) as i32;
            }
            b_min = if b_min == -1.0 { major as f64 } else { jmin(major as f64, b_min) };
        }
    }
    for s in &mut allele[start_index..=end_index] {
        s.a_fixed_ap = a_jcn;
        s.b_undisrupted_ap = b_min;
        s.set_valid(true);
    }
    true
}

pub fn has_valid_allele_jcn_data(breakend_index: i32, segments: Option<&Vec<SegmentJcn>>) -> bool {
    let Some(s) = segments else { return false };
    if (s.len() as i32) < breakend_index {
        return false;
    }
    // Java: List.get throws for an index equal to the size
    s[breakend_index as usize].is_valid()
}

/// ChainJcnLimits.calcJcnUncertainty
pub fn calc_jcn_uncertainty(d1: &JcnCalcData, d2: &JcnCalcData) -> JcnCalcData {
    if !d1.valid || !d2.valid {
        return JcnCalcData { jcn_estimate: 0.0, jcn_uncertainty: 0.0, valid: false };
    }
    if d1.jcn_uncertainty == 0.0 {
        return JcnCalcData { jcn_estimate: d1.jcn_estimate, jcn_uncertainty: 0.0, valid: true };
    } else if d2.jcn_uncertainty == 0.0 {
        return JcnCalcData { jcn_estimate: d2.jcn_estimate, jcn_uncertainty: 0.0, valid: true };
    }
    let u1 = 1.0 / d1.jcn_uncertainty.powi(2);
    let u2 = 1.0 / d2.jcn_uncertainty.powi(2);
    let sum_u = u1 + u2;
    let sum_ou = d1.jcn_estimate * u1 + d2.jcn_estimate * u2;
    let est = sum_ou / sum_u;
    let mut adj = u1 * jmax(d1.jcn_estimate - est, d1.jcn_uncertainty / 2.0).powi(2);
    adj += u2 * jmax(d2.jcn_estimate - est, d2.jcn_uncertainty / 2.0).powi(2);
    let est_u = (2.0 * adj / sum_u).sqrt();
    JcnCalcData { jcn_estimate: est, jcn_uncertainty: est_u, valid: true }
}

pub const UNCERTAINTY_SCALE_FACTOR: f64 = std::f64::consts::SQRT_2;

pub fn jcn_match_for_splits(jcn_lower: f64, unc_lower: f64, jcn: f64, unc: f64) -> bool {
    copy_numbers_equal(jcn_lower * 2.0, jcn) || jcn_overlap(jcn_lower * 2.0, unc_lower * UNCERTAINTY_SCALE_FACTOR, jcn, unc)
}

pub fn jcn_match_svs(var1: &SvVarData, var2: &SvVarData) -> bool {
    jcn_match(var1.jcn(), var1.jcn_uncertainty(), var2.jcn(), var2.jcn_uncertainty())
}

pub fn jcn_match_breakends(b1: &SvBreakend, b2: &SvBreakend) -> bool {
    jcn_match(b1.jcn(), b1.jcn_uncertainty(), b2.jcn(), b2.jcn_uncertainty())
}

pub fn jcn_match(jcn1: f64, unc1: f64, jcn2: f64, unc2: f64) -> bool {
    copy_numbers_equal(jcn1, jcn2) || jcn_overlap(jcn1, unc1, jcn2, unc2)
}

pub fn jcn_overlap(jcn1: f64, unc1: f64, jcn2: f64, unc2: f64) -> bool {
    !(jcn1 + unc1 < jcn2 - unc2 || jcn2 + unc2 < jcn1 - unc1)
}

pub fn jcn_exceeds_major_allele_jcn(breakend: &SvBreakend) -> bool {
    let adjacent = breakend.major_allele_jcn(breakend.reverse_orient());
    breakend.jcn() > adjacent && !copy_numbers_equal(breakend.jcn(), adjacent)
}
