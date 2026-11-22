//! AssemblyLink: two assemblies joined as a split variant, an indel or a facing pair.

use crate::assembly::assembly_utils::find_matching_assembly;
use crate::assembly::ptr::P;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::{form_sv_type, LinkType, StructuralVariantType};

pub type LinkP = P<AssemblyLink>;

pub struct AssemblyLink {
    first: AsmP,
    second: AsmP,
    type_: LinkType,
    inserted_bases: String,
    overlap_bases: String,
    insert_site: bool,
}

impl AssemblyLink {
    pub fn new(first: AsmP, second: AsmP, type_: LinkType, inserted_bases: &str, overlap_bases: &str) -> LinkP {
        P::new(AssemblyLink { first, second, type_, inserted_bases: inserted_bases.to_string(), overlap_bases: overlap_bases.to_string(), insert_site: false })
    }
    pub fn type_(&self) -> LinkType {
        self.type_
    }
    pub fn first(&self) -> AsmP {
        self.first
    }
    pub fn second(&self) -> AsmP {
        self.second
    }
    pub fn mark_insert_site(&mut self) {
        self.insert_site = true;
    }
    pub fn is_insert_site(&self) -> bool {
        self.insert_site
    }
    pub fn other_assembly(&self, assembly: AsmP) -> AsmP {
        if self.first == assembly { self.second } else { self.first }
    }
    pub fn inserted_bases(&self) -> &str {
        &self.inserted_bases
    }
    pub fn overlap_bases(&self) -> &str {
        &self.overlap_bases
    }
    pub fn sv_type(&self) -> StructuralVariantType {
        let (j1, j2) = (self.first.junction(), self.second.junction());
        form_sv_type(j1.chromosome, j2.chromosome, j1.position, j2.position, j1.orient, j2.orient, !self.inserted_bases.is_empty())
    }
    pub fn matches_link(&self, other: &AssemblyLink) -> bool {
        self.matches(other.first, other.second)
    }
    pub fn matches(&self, assembly1: AsmP, assembly2: AsmP) -> bool {
        (self.first == assembly1 && self.second == assembly2) || (self.first == assembly2 && self.second == assembly1)
    }
    pub fn has_assembly(&self, assembly: AsmP) -> bool {
        self.first == assembly || self.second == assembly
    }
    pub fn find_matched_assembly(&self, assembly: AsmP) -> Option<AsmP> {
        find_matching_assembly(&[self.first, self.second], assembly, true)
    }
    pub fn length(&self) -> i32 {
        let (j1, j2) = (self.first.junction(), self.second.junction());
        if self.type_ == LinkType::Facing {
            return (j1.position - j2.position).abs() + 1;
        }
        if j1.chromosome == j2.chromosome { (j1.position - j2.position).abs() } else { 0 }
    }
    pub fn swap_assemblies(existing_link: LinkP, original_assembly: AsmP, new_assembly: AsmP) -> LinkP {
        let (first, second) = if existing_link.first == original_assembly { (new_assembly, existing_link.second) } else { (existing_link.first, new_assembly) };
        AssemblyLink::new(first, second, existing_link.type_, &existing_link.inserted_bases, &existing_link.overlap_bases)
    }
}

impl std::fmt::Display for AssemblyLink {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}: {} - {} len({}) extras(overlap={} insert={})",
            self.type_.name(),
            self.first.junction().coords(),
            self.second.junction().coords(),
            self.length(),
            self.overlap_bases.len(),
            self.inserted_bases.len()
        )
    }
}
