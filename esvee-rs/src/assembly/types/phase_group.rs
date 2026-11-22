//! PhaseGroup: assemblies that share reads or proximity and are linked together.

use crate::assembly::alignment::assembly_alignment::AssemblyAlignment;
use crate::assembly::phase::phase_set_merger;
use crate::assembly::ptr::P;
use crate::assembly::types::assembly_link::LinkP;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::phase_set::PhaseSetP;
use crate::assembly::types::AssemblyOutcome;

pub type PhaseGroupP = P<PhaseGroup>;

pub struct PhaseGroup {
    id: i32,
    assemblies: Vec<AsmP>,
    derived_assemblies: Vec<AsmP>,
    phase_sets: Vec<PhaseSetP>,
    secondary_split_links: Vec<LinkP>,
    next_phase_set_id: i32,
}

impl PhaseGroup {
    pub fn new(first: AsmP, second: Option<AsmP>) -> PhaseGroupP {
        let pg = P::new(PhaseGroup { id: -1, assemblies: vec![first], derived_assemblies: Vec::new(), phase_sets: Vec::new(), secondary_split_links: Vec::new(), next_phase_set_id: 0 });
        first.m().set_phase_group(Some(pg));
        if let Some(second) = second {
            pg.m().assemblies.push(second);
            second.m().set_phase_group(Some(pg));
        }
        pg
    }
    pub fn set_id(&mut self, id: i32) {
        self.id = id;
    }
    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn assemblies(&self) -> &Vec<AsmP> {
        &self.assemblies
    }
    pub fn assemblies_mut(&mut self) -> &mut Vec<AsmP> {
        &mut self.assemblies
    }
    pub fn derived_assemblies(&self) -> &Vec<AsmP> {
        &self.derived_assemblies
    }
    pub fn derived_assemblies_mut(&mut self) -> &mut Vec<AsmP> {
        &mut self.derived_assemblies
    }
    pub fn phase_sets(&self) -> &Vec<PhaseSetP> {
        &self.phase_sets
    }
    pub fn phase_sets_mut(&mut self) -> &mut Vec<PhaseSetP> {
        &mut self.phase_sets
    }
    pub fn secondary_links(&self) -> &Vec<LinkP> {
        &self.secondary_split_links
    }
    pub fn secondary_links_mut(&mut self) -> &mut Vec<LinkP> {
        &mut self.secondary_split_links
    }
    pub fn find_phase_set(&self, assembly: AsmP) -> Option<PhaseSetP> {
        for phase_set in &self.phase_sets {
            if phase_set.has_assembly(assembly) {
                return Some(*phase_set);
            }
            if assembly.outcome() == AssemblyOutcome::Secondary && phase_set.secondary_links().iter().any(|x| x.has_assembly(assembly)) {
                return Some(*phase_set);
            }
        }
        None
    }
    pub fn find_secondary_split_links(&self, assembly: AsmP) -> Vec<LinkP> {
        self.secondary_split_links.iter().filter(|x| x.has_assembly(assembly)).copied().collect()
    }
    pub fn assembly_count(&self) -> usize {
        self.assemblies.len()
    }

    pub fn transfer_assemblies(this: PhaseGroupP, other: PhaseGroupP) {
        for assembly in other.assemblies.iter() {
            if this.assemblies.contains(assembly) {
                crate::error!("assembly({}) transferred from pg({})but already in phase group({})", **assembly, *other, *this);
                std::process::exit(1);
            }
            this.m().assemblies.push(*assembly);
            assembly.m().set_phase_group(Some(this));
        }
    }

    pub fn add_assembly(this: PhaseGroupP, assembly: AsmP) {
        this.m().assemblies.push(assembly);
        if assembly.phase_group().is_some() {
            crate::error!("assembly({}) adding to additional phase group", *assembly);
            std::process::exit(1);
        }
        assembly.m().set_phase_group(Some(this));
    }

    pub fn add_derived_assembly(this: PhaseGroupP, assembly: AsmP) {
        this.m().derived_assemblies.push(assembly);
        Self::add_assembly(this, assembly);
    }

    pub fn finalise_phase_set_alignments(this: PhaseGroupP) {
        for phase_set in this.phase_sets.clone() {
            let id = this.m().next_phase_set_id();
            phase_set.m().set_id(id);
            if phase_set.is_short_local_ref_link() {
                continue;
            }
            let assembly_alignment = AssemblyAlignment::from_phase_set(phase_set);
            phase_set.m().set_assembly_alignment(Some(assembly_alignment));
        }
        phase_set_merger::merge_phase_sets(this.m().phase_sets_mut());
    }

    pub fn next_phase_set_id(&mut self) -> i32 {
        let id = self.next_phase_set_id;
        self.next_phase_set_id += 1;
        id
    }
}

impl std::fmt::Display for PhaseGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "id({}) assemblies({})", self.id, self.assemblies.len())
    }
}
