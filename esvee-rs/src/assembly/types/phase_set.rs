//! PhaseSet: a chain of linked assemblies.

use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::constants::{MIN_VARIANT_LENGTH, MIN_VARIANT_LENGTH_LOWER};
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::P;
use crate::assembly::types::assembly_link::LinkP;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::{AssemblyOutcome, LinkType};

pub type PhaseSetP = P<PhaseSet>;

pub struct PhaseSet {
    id: i32,
    assembly_links: Vec<LinkP>,
    secondary_links: Vec<LinkP>,
    assemblies: Vec<AsmP>,
    assembly_alignment: Option<AsmAlignP>,
    merged_phase_sets: Vec<PhaseSetP>,
    merged_phase_set_id: Option<i32>,
}

impl PhaseSet {
    pub fn new(link: LinkP) -> PhaseSetP {
        let mut ps = PhaseSet {
            id: -1,
            assembly_links: Vec::new(),
            secondary_links: Vec::new(),
            assemblies: Vec::new(),
            assembly_alignment: None,
            merged_phase_sets: Vec::new(),
            merged_phase_set_id: None,
        };
        ps.add_assembly_link(link, 0);
        P::new(ps)
    }
    pub fn set_id(&mut self, id: i32) {
        self.id = id;
    }
    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn add_assembly_link_start(&mut self, link: LinkP) {
        self.add_assembly_link(link, 0);
    }
    pub fn add_assembly_link_end(&mut self, link: LinkP) {
        let n = self.assembly_links.len();
        self.add_assembly_link(link, n);
    }
    pub fn add_secondary_link(&mut self, link: LinkP) {
        self.secondary_links.push(link);
    }
    fn add_assembly_link(&mut self, link: LinkP, index: usize) {
        self.assembly_links.insert(index, link);
        if self.assemblies.is_empty() {
            self.assemblies.push(link.first());
            self.assemblies.push(link.second());
            return;
        }
        let assembly_index = if index == 0 { 0 } else { self.assemblies.len() };
        if !self.assemblies.contains(&link.first()) {
            self.assemblies.insert(assembly_index, link.first());
        }
        if !self.assemblies.contains(&link.second()) {
            self.assemblies.insert(assembly_index, link.second());
        }
    }
    pub fn assemblies(&self) -> &Vec<AsmP> {
        &self.assemblies
    }
    pub fn all_assemblies(&self) -> Vec<AsmP> {
        if self.secondary_links.is_empty() {
            return self.assemblies.clone();
        }
        let mut assemblies = self.assemblies.clone();
        assemblies.extend(self.secondary_assemblies());
        assemblies
    }
    pub fn secondary_assemblies(&self) -> Vec<AsmP> {
        let mut assemblies: Vec<AsmP> = Vec::new();
        for link in &self.secondary_links {
            let first = link.first();
            if !assemblies.contains(&first) && !self.assemblies.contains(&first) {
                assemblies.push(first);
            }
            let second = link.second();
            if !assemblies.contains(&second) && !self.assemblies.contains(&second) {
                assemblies.push(second);
            }
        }
        assemblies
    }
    pub fn assembly_links(&self) -> &Vec<LinkP> {
        &self.assembly_links
    }
    pub fn secondary_links(&self) -> &Vec<LinkP> {
        &self.secondary_links
    }
    pub fn has_valid_assembly_alignment(&self) -> bool {
        self.assembly_alignment.is_some_and(|a| a.is_valid())
    }
    pub fn assembly_alignment(&self) -> Option<AsmAlignP> {
        self.assembly_alignment
    }
    pub fn set_assembly_alignment(&mut self, assembly_alignment: Option<AsmAlignP>) {
        self.assembly_alignment = assembly_alignment;
    }
    pub fn merged_phase_sets(&self) -> &Vec<PhaseSetP> {
        &self.merged_phase_sets
    }
    pub fn merge_phase_set(&mut self, phase_set: PhaseSetP) {
        self.merged_phase_sets.push(phase_set);
    }
    pub fn merged(&self) -> bool {
        self.merged_phase_set_id.is_some()
    }
    pub fn merged_phase_set_id(&self) -> i32 {
        self.merged_phase_set_id.unwrap_or(-1)
    }
    pub fn set_merged_phase_set_id(&mut self, phase_set_id: i32) {
        self.merged_phase_set_id = Some(phase_set_id);
    }
    pub fn has_assembly(&self, assembly: AsmP) -> bool {
        self.assemblies.contains(&assembly)
    }
    pub fn find_assembly_links(&self, assembly: AsmP) -> Vec<LinkP> {
        self.assembly_links.iter().filter(|x| x.has_assembly(assembly)).copied().collect()
    }
    pub fn find_split_link(&self, assembly: AsmP) -> Option<LinkP> {
        self.assembly_links.iter().filter(|x| x.has_assembly(assembly)).find(|x| x.type_() == LinkType::Split).copied()
    }
    pub fn is_secondary_line_link(&self) -> bool {
        self.assembly_links.len() == 1 && self.assemblies.iter().any(|x| x.outcome() == AssemblyOutcome::Secondary)
    }
    pub fn is_short_local_ref_link(&self) -> bool {
        if self.assembly_links.len() != 1 || self.assemblies[0].outcome() != AssemblyOutcome::LocalIndel {
            return false;
        }
        let is_saga_matched = self.assemblies.iter().any(|x| x.is_saga_matched());
        let threshold = if is_saga_matched { MIN_VARIANT_LENGTH_LOWER } else { MIN_VARIANT_LENGTH };
        self.assembly_links[0].length() < threshold
    }
    pub fn has_facing_links(&self) -> bool {
        self.assembly_links.iter().any(|x| x.type_() == LinkType::Facing)
    }
    pub fn assemblies_face_in_phase_set(&self, assembly1: AsmP, assembly2: AsmP) -> bool {
        let mut assembly_orientation = Orientation::Forward;
        let mut assembly_index1 = -1;
        let mut assembly_orientation1: Option<Orientation> = None;
        let mut assembly_index2 = -1;
        let mut assembly_orientation2: Option<Orientation> = None;
        for (i, a) in self.assemblies.iter().enumerate() {
            if *a == assembly1 {
                assembly_index1 = i as i32;
                assembly_orientation1 = Some(assembly_orientation);
                if assembly_orientation2.is_some() {
                    break;
                }
            } else if *a == assembly2 {
                assembly_index2 = i as i32;
                assembly_orientation2 = Some(assembly_orientation);
                if assembly_orientation1.is_some() {
                    break;
                }
            }
            assembly_orientation = assembly_orientation.opposite();
        }
        let (Some(o1), Some(o2)) = (assembly_orientation1, assembly_orientation2) else { return false };
        if o1 == o2 {
            return false;
        }
        if assembly_index1 < assembly_index2 { o1.is_forward() } else { o2.is_forward() }
    }
}

impl std::fmt::Display for PhaseSet {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.assembly_links.len() == 1 {
            return write!(f, "id({}) link({})", self.id, *self.assembly_links[0]);
        }
        write!(
            f,
            "id({}) links({} facing={}, seconds={}) first({})",
            self.id,
            self.assembly_links.len(),
            self.assembly_links.iter().filter(|x| x.type_() == LinkType::Facing).count(),
            self.secondary_links.len(),
            *self.assembly_links[0]
        )
    }
}
