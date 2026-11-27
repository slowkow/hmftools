//! LineChainer: chains of a LINE cluster's source breakends (its cohort output is left out).

use std::collections::HashSet;

use crate::annotators::LINE_ELEMENT_PROXIMITY_DISTANCE;
use crate::chaining::utils::{get_min_templated_insertion_length, is_possible_link};
use crate::types::*;

#[derive(Default)]
pub struct LineChainer {
    cluster_id: i32,
    sample_id: String,
    sv_list: Vec<Sv>,
    assembled_links: Vec<Pair>,
    chains: Vec<Chain>,
    source_breakends: Vec<Be>,
    source_chromosomes: HashSet<String>,
}

impl LineChainer {
    pub fn new() -> Self {
        LineChainer { cluster_id: -1, ..Default::default() }
    }
    pub fn chains(&self) -> &[Chain] {
        &self.chains
    }
    pub fn clear(&mut self) {
        self.cluster_id = -1;
        self.sv_list.clear();
        self.chains.clear();
        self.assembled_links.clear();
        self.source_breakends.clear();
        self.source_chromosomes.clear();
    }
    pub fn initialise(&mut self, sample_id: &str, cluster: &Cluster) {
        self.clear();
        self.cluster_id = cluster.id();
        self.sample_id = sample_id.to_string();
        self.sv_list = cluster.svs();
        self.assembled_links = cluster.assembly_linked_pairs();
    }

    pub fn form_chains(&mut self) {
        let mut sgl_mapped: Vec<Sv> = Vec::new();
        for var in self.sv_list.clone() {
            for s in [SE_START, SE_END] {
                let Some(breakend) = var.breakend_se(s) else {
                    if !var.sgl_mappings().is_empty() {
                        sgl_mapped.push(var.clone());
                    }
                    continue;
                };
                if !var.is_line_element(is_start(s)) {
                    continue;
                }
                self.source_chromosomes.insert(breakend.chromosome().to_string());
                if self.assembled_links.iter().any(|x| x.has_breakend(&breakend)) {
                    continue;
                }
                self.source_breakends.push(breakend);
            }
        }
        self.add_sgl_mappings(&sgl_mapped);
        for pair in self.assembled_links.clone() {
            if self.try_add_link(&pair) {
                continue;
            }
            self.add_new_chain(&pair);
        }
        if !self.source_breakends.is_empty() {
            self.add_source_breakends();
        }
    }

    fn add_new_chain(&mut self, pair: &Pair) {
        let chain = Chain::create(self.chains.len() as i32);
        for s in [SE_START, SE_END] {
            let b = pair.breakend_se(s);
            if b.sv().is_sgl_breakend() && !b.uses_start() {
                b.sv().set_sgl_end_breakend(Some(b.clone()));
            }
        }
        chain.add_link(pair, true);
        hmf_common::debug!("new chain({}) with pair({})", chain.id(), pair);
        self.chains.push(chain);
    }

    fn try_add_link(&self, pair: &Pair) -> bool {
        for chain in &self.chains {
            for s in [SE_START, SE_END] {
                let Some(cb) = chain.open_breakend(is_start(s)) else { continue };
                if !self.source_chromosomes.contains(cb.chromosome()) {
                    continue;
                }
                if chain.can_add_linked_pair(pair, is_start(s), false) {
                    chain.add_link(pair, is_start(s));
                    return true;
                }
            }
        }
        false
    }

    fn add_source_breakends(&mut self) {
        loop {
            let mut shortest: Option<Pair> = None;
            for i in 0..self.source_breakends.len() {
                let b = self.source_breakends[i].clone();
                for j in i + 1..self.source_breakends.len() {
                    let b2 = self.source_breakends[j].clone();
                    let Some(p) = self.try_form_linked_pair(&b, &b2) else { continue };
                    if shortest.as_ref().is_none_or(|s| p.position_distance() < s.position_distance()) {
                        shortest = Some(p);
                    }
                }
                for chain in self.chains.clone() {
                    for s in [SE_START, SE_END] {
                        let Some(cb) = chain.open_breakend(is_start(s)) else { continue };
                        if !self.source_chromosomes.contains(cb.chromosome()) {
                            continue;
                        }
                        let Some(p) = self.try_form_linked_pair(&b, &cb) else { continue };
                        if shortest.as_ref().is_none_or(|s| p.position_distance() < s.position_distance()) {
                            shortest = Some(p);
                        }
                    }
                }
            }
            let Some(sp) = shortest else { return };
            if !self.try_add_link(&sp) {
                self.add_new_chain(&sp);
            }
            remove_item(&mut self.source_breakends, &sp.first_breakend());
            remove_item(&mut self.source_breakends, &sp.second_breakend());
            if self.source_breakends.is_empty() {
                return;
            }
        }
    }

    fn try_form_linked_pair(&self, b1: &Be, b2: &Be) -> Option<Pair> {
        let min_ti = get_min_templated_insertion_length(b2, b1);
        if !is_possible_link(b1.chromosome(), b1.position(), b1.orientation(), b2.chromosome(), b2.position(), b2.orientation(), min_ti) {
            return None;
        }
        if (b1.position() - b2.position()).abs() > LINE_ELEMENT_PROXIMITY_DISTANCE * 2 {
            return None;
        }
        let o1 = self.other_breakend(b1)?;
        let o2 = self.other_breakend(b2)?;
        let d1 = o1.db_link();
        if d1.is_none() || d1 != o2.db_link() {
            return None;
        }
        Some(Pair::from_breakends(b1, b2))
    }

    fn other_breakend(&self, breakend: &Be) -> Option<Be> {
        for chain in &self.chains {
            for s in [SE_START, SE_END] {
                if chain.open_breakend(is_start(s)).as_ref() == Some(breakend) {
                    return chain.open_breakend(!is_start(s));
                }
            }
        }
        breakend.other_breakend()
    }

    fn add_sgl_mappings(&mut self, sgl_mapped: &[Sv]) {
        let mut candidates: Vec<Sv> = Vec::new();
        if self.source_breakends.is_empty() {
            for i in 0..sgl_mapped.len() {
                let sgl1 = &sgl_mapped[i];
                if contains(&candidates, sgl1) {
                    continue;
                }
                for sgl2 in sgl_mapped.iter().skip(i + 1) {
                    if contains(&candidates, sgl2) {
                        continue;
                    }
                    let m2s = sgl2.sgl_mappings();
                    for m1 in sgl1.sgl_mappings() {
                        let m2 = m2s.iter().find(|x| {
                            is_possible_link(
                                &x.chromosome,
                                x.position,
                                x.orientation,
                                &m1.chromosome,
                                m1.position,
                                m1.orientation,
                                MIN_TEMPLATED_INSERTION_LENGTH,
                            )
                        });
                        if let Some(m2) = m2 {
                            self.source_breakends.push(Be::new(SvBreakend::from_mapping(sgl1, &m1)));
                            self.source_breakends.push(Be::new(SvBreakend::from_mapping(sgl2, m2)));
                            candidates.push(sgl1.clone());
                            candidates.push(sgl2.clone());
                            sgl1.add_line_element(LineElementType::SUSPECT, false);
                            sgl2.add_line_element(LineElementType::SUSPECT, false);
                            break;
                        }
                    }
                }
            }
        } else {
            for sgl in sgl_mapped {
                let mut linked = false;
                for bi in 0..self.source_breakends.len() {
                    let b = self.source_breakends[bi].clone();
                    let min_distance = (b.anchor_distance() - b.homology().len() as i32).max(MIN_TEMPLATED_INSERTION_LENGTH);
                    for m in sgl.sgl_mappings() {
                        if is_possible_link(&m.chromosome, m.position, m.orientation, b.chromosome(), b.position(), b.orientation(), min_distance) {
                            self.source_breakends.push(Be::new(SvBreakend::from_mapping(sgl, &m)));
                            sgl.add_line_element(LineElementType::SUSPECT, false);
                            linked = true;
                            break;
                        }
                    }
                    if linked {
                        break;
                    }
                }
            }
        }
    }
}
