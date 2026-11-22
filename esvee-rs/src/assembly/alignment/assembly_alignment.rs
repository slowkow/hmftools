//! AssemblyAlignment: the full sequence formed from a single assembly or from the linked assemblies of a
//! phase set, which is what gets aligned to the reference.

use std::collections::HashMap;

use crate::assembly::alignment::breakend::BreakendP;
use crate::assembly::assembly_utils::{extract_insert_sequence, has_unset_bases, reverse_complement_str};
use crate::assembly::constants::ASSEMBLY_INFO_DELIM;
use crate::assembly::globals::assembly_build_debug;
use crate::assembly::orientation::Orientation;
use crate::assembly::phase::assembly_linker::is_assembly_indel_link;
use crate::assembly::ptr::P;
use crate::assembly::saga::SagaSequenceMatch;
use crate::assembly::sam::{cigar_elements_to_str, CigarElement, CigarOp};
use crate::assembly::types::assembly_link::LinkP;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::phase_set::PhaseSetP;
use crate::assembly::types::support_read::SupP;
use crate::assembly::types::{AssemblyOutcome, LinkType};

pub type AsmAlignP = P<AssemblyAlignment>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum ReadLinkType {
    Left,
    Right,
    None,
}

pub struct AssemblyAlignment {
    id: i32,
    /// the single assembly, when not built from a phase set
    single_assembly: Option<AsmP>,
    phase_set: Option<PhaseSetP>,
    full_sequence: Option<String>,
    full_sequence_length: i32,
    sequence_cigar: String,
    breakends: Vec<BreakendP>,
    link_indices: Vec<i32>,
    fragment_reads_map: HashMap<String, Vec<SupP>>,
    saga_match: Option<SagaSequenceMatch>,
}

impl AssemblyAlignment {
    pub fn from_assembly(assembly: AsmP) -> AsmAlignP {
        Self::build(Some(assembly), None)
    }
    pub fn from_phase_set(phase_set: PhaseSetP) -> AsmAlignP {
        Self::build(None, Some(phase_set))
    }

    fn build(assembly: Option<AsmP>, phase_set: Option<PhaseSetP>) -> AsmAlignP {
        let mut aa = AssemblyAlignment {
            id: -1,
            single_assembly: assembly,
            phase_set,
            full_sequence: None,
            full_sequence_length: 0,
            sequence_cigar: String::new(),
            breakends: Vec::new(),
            link_indices: Vec::new(),
            fragment_reads_map: HashMap::new(),
            saga_match: None,
        };
        let info = aa.info();
        for x in aa.assemblies() {
            x.m().set_assembly_alignment_info(Some(info.clone()));
        }
        if !aa.assemblies().iter().any(|x| has_unset_bases(x)) {
            aa.full_sequence = aa.build_sequence_data();
            aa.fragment_reads_map.clear();
        }
        aa.full_sequence_length = aa.full_sequence.as_ref().map(|s| s.len() as i32).unwrap_or(0);
        P::new(aa)
    }

    pub fn set_id(&mut self, id: i32) {
        self.id = id;
    }
    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn is_valid(&self) -> bool {
        self.full_sequence.is_some() && self.full_sequence_length < 100000
    }
    /// the phase set's (live) list of assemblies, or the single assembly
    pub fn assemblies(&self) -> Vec<AsmP> {
        match (&self.single_assembly, &self.phase_set) {
            (Some(a), _) => vec![*a],
            (None, Some(ps)) => ps.assemblies().clone(),
            _ => Vec::new(),
        }
    }
    pub fn phase_set(&self) -> Option<PhaseSetP> {
        self.phase_set
    }
    pub fn is_merged(&self) -> bool {
        self.phase_set.is_some_and(|ps| ps.merged())
    }
    pub fn breakends(&self) -> &Vec<BreakendP> {
        &self.breakends
    }
    pub fn breakends_mut(&mut self) -> &mut Vec<BreakendP> {
        &mut self.breakends
    }
    pub fn add_breakend(&mut self, breakend: BreakendP) {
        self.breakends.push(breakend);
    }
    /// null in Java when an assembly has unset bases or the sequence could not be built
    pub fn full_sequence(&self) -> Option<&str> {
        self.full_sequence.as_deref()
    }
    pub fn full_sequence_length(&self) -> i32 {
        self.full_sequence_length
    }
    pub fn assembly_cigar(&self) -> &str {
        &self.sequence_cigar
    }
    pub fn link_indices(&self) -> &Vec<i32> {
        &self.link_indices
    }
    pub fn assembly_ids(&self) -> String {
        self.assemblies().iter().map(|x| x.id().to_string()).collect::<Vec<_>>().join(";")
    }
    pub fn info(&self) -> String {
        self.info_limit(-1)
    }
    pub fn info_limit(&self, limit: i32) -> String {
        let assemblies = self.assemblies();
        if assemblies.len() == 1 {
            return assemblies[0].junction().coords_typed_ext(true);
        }
        let mut parts: Vec<String> = Vec::new();
        for (i, a) in assemblies.iter().enumerate() {
            if limit > 0 && i as i32 >= limit {
                break;
            }
            parts.push(a.junction().coords_typed_ext(true));
        }
        parts.join(ASSEMBLY_INFO_DELIM)
    }

    fn build_single_assembly_sequence_data(&mut self) -> Option<String> {
        let assembly = self.assemblies()[0];
        let extension_length = assembly.extension_length();
        let ref_base_length = assembly.ref_base_length();
        if assembly.is_forward_junction() {
            self.set_assembly_read_indices(assembly, false, ref_base_length - 1, ReadLinkType::None);
            self.sequence_cigar = format!("{ref_base_length}M{extension_length}S");
            self.link_indices.push(ref_base_length);
        } else {
            self.set_assembly_read_indices(assembly, false, extension_length, ReadLinkType::None);
            self.sequence_cigar = format!("{extension_length}S{ref_base_length}M");
            self.link_indices.push(extension_length);
        }
        Some(assembly.form_full_sequence())
    }

    fn build_sequence_data(&mut self) -> Option<String> {
        let mut assemblies = self.assemblies();
        if assemblies.len() == 1 {
            return self.build_single_assembly_sequence_data();
        }
        let phase_set = self.phase_set.expect("phase set");
        let mut assembly_links: Vec<LinkP> = phase_set.assembly_links().clone();
        let assembly_count = assemblies.len();
        let link_count = assembly_links.len();
        let mut has_facing_at_start = false;
        let mut has_facing_at_end = false;
        if phase_set.has_facing_links() {
            has_facing_at_start = assembly_links[0].type_() == LinkType::Facing;
            has_facing_at_end = assembly_links[link_count - 1].type_() == LinkType::Facing;
        }
        let has_outer_extensions = has_facing_at_start && has_facing_at_end;
        let reverse_links = has_facing_at_start && !has_facing_at_end;
        if reverse_links {
            has_facing_at_end = true;
            assemblies.reverse();
            assembly_links.reverse();
        }
        let start_reversed = if has_outer_extensions { false } else { assemblies[0].junction().is_reverse() };
        let mut full_sequence = String::new();
        let mut last_added_reversed = false;
        let mut current_seq_length: i32 = 0;
        let mut sequence_cigar: Vec<CigarElement> = Vec::new();

        for i in 0..assembly_count - 1 {
            let link = assembly_links[i];
            let assembly = assemblies[i];
            let is_facing_link = link.type_() == LinkType::Facing;
            let assembly_reversed;
            if i == 0 {
                assembly_reversed = start_reversed;
                if is_facing_link {
                    let assembly_extension_bases = assembly.form_junction_sequence();
                    let assembly_extension_base_length = assembly_extension_bases.len() as i32;
                    full_sequence.push_str(&assembly_extension_bases);
                    self.set_assembly_read_indices(assembly, assembly_reversed, assembly_extension_base_length, ReadLinkType::None);
                    log_build_info(assembly, current_seq_length, assembly_extension_base_length, assembly_reversed, "outer-ext-bases");
                    current_seq_length = assembly_extension_base_length;
                    build_sequence_cigar(&mut sequence_cigar, CigarOp::S, assembly_extension_base_length);
                    self.link_indices.push(assembly_extension_base_length);
                } else {
                    let assembly_ref_bases = if start_reversed { reverse_complement_str(&assembly.form_ref_base_sequence()) } else { assembly.form_ref_base_sequence() };
                    let assembly_ref_base_length = assembly_ref_bases.len() as i32;
                    full_sequence.push_str(&assembly_ref_bases);
                    self.set_assembly_read_indices(assembly, assembly_reversed, assembly_ref_base_length - 1, ReadLinkType::Left);
                    current_seq_length += assembly_ref_base_length;
                    build_sequence_cigar(&mut sequence_cigar, CigarOp::M, assembly_ref_base_length);
                    log_build_info(assembly, current_seq_length, assembly_ref_base_length, assembly_reversed, "ref-bases");
                }
            } else {
                if is_facing_link {
                    let next_assembly = assemblies[i + 1];
                    self.set_assembly_read_indices(next_assembly, last_added_reversed, current_seq_length, ReadLinkType::Left);
                    if i == assembly_count - 2 {
                        add_final_facing_extension_bases(next_assembly, &mut full_sequence, &mut sequence_cigar, current_seq_length, last_added_reversed, &assembly_links);
                        current_seq_length += next_assembly.extension_length();
                    }
                    continue;
                }
                assembly_reversed = last_added_reversed;
            }
            let overlap_length = link.overlap_bases().len() as i32;
            let inserted_base_length = link.inserted_bases().len() as i32;
            let next_assembly = assemblies[i + 1];
            let next_reversed = next_assembly.is_forward_junction() && link.type_() != LinkType::Facing;
            if inserted_base_length > 0 {
                let inserted_bases = if link.is_insert_site() {
                    link.inserted_bases().to_string()
                } else {
                    extract_insert_sequence(&assembly, assembly_reversed, &next_assembly, next_reversed, inserted_base_length)
                };
                full_sequence.push_str(&inserted_bases);
                current_seq_length += inserted_base_length;
                build_sequence_cigar(&mut sequence_cigar, CigarOp::I, inserted_base_length);
            }
            let mut next_assembly_ref_bases = if next_reversed { reverse_complement_str(&next_assembly.form_ref_base_sequence()) } else { next_assembly.form_ref_base_sequence() };
            if overlap_length >= next_assembly_ref_bases.len() as i32 {
                return None;
            }
            if overlap_length > 0 {
                next_assembly_ref_bases = next_assembly_ref_bases[overlap_length as usize..].to_string();
            }
            full_sequence.push_str(&next_assembly_ref_bases);
            self.link_indices.push(current_seq_length);
            let next_len = next_assembly_ref_bases.len() as i32;
            let next_assembly_junction_index = if link.type_() != LinkType::Facing { current_seq_length } else { current_seq_length + next_len };
            self.set_assembly_read_indices(next_assembly, next_reversed, next_assembly_junction_index, ReadLinkType::Right);
            log_build_info(next_assembly, current_seq_length, next_len, assembly_reversed, "ref-bases");
            if overlap_length > 0 {
                build_sequence_cigar(&mut sequence_cigar, CigarOp::D, overlap_length);
            }
            current_seq_length += next_len;
            build_sequence_cigar(&mut sequence_cigar, CigarOp::M, next_len);
            if has_facing_at_end && i == assembly_count - 2 {
                add_final_facing_extension_bases(next_assembly, &mut full_sequence, &mut sequence_cigar, current_seq_length, assembly_reversed, &assembly_links);
                // as in Java, assigned rather than added
                current_seq_length = next_assembly.extension_length();
            }
            last_added_reversed = next_reversed;
        }
        let _ = current_seq_length;
        self.sequence_cigar = cigar_elements_to_str(&sequence_cigar);
        Some(full_sequence)
    }

    fn set_assembly_read_indices(&mut self, assembly: AsmP, is_reversed: bool, assembly_junction_seq_index: i32, link_type: ReadLinkType) {
        let current_seq_end_index = assembly_junction_seq_index;
        for read in assembly.support().clone() {
            let junc_read_start_distance = read.junction_read_start_distance();
            let mut full_seq_orientation = if is_reversed { Orientation::Reverse } else { Orientation::Forward };
            let mut full_seq_index = current_seq_end_index;
            if !is_reversed {
                full_seq_index -= junc_read_start_distance;
            } else {
                full_seq_index += junc_read_start_distance - read.base_length() + 1;
            }
            if let Some(matched_read) = self.find_matching_support_read(&read) {
                full_seq_orientation = matched_read.full_assembly_orientation();
                let use_matched_read_info;
                let mut update_matched_read_info = true;
                if read.trim_count_start() != matched_read.trim_count_start() || read.trim_count_end() != matched_read.trim_count_end() {
                    use_matched_read_info = false;
                    update_matched_read_info = false;
                } else if read.type_() == matched_read.type_() && !read.type_().is_split_support() {
                    use_matched_read_info = junc_read_start_distance > matched_read.junction_read_start_distance();
                } else {
                    use_matched_read_info = true;
                }
                if use_matched_read_info {
                    full_seq_index = matched_read.full_assembly_index_start();
                } else if update_matched_read_info {
                    matched_read.m().set_full_assembly_info(full_seq_index, full_seq_orientation);
                }
            }
            read.m().set_full_assembly_info(full_seq_index, full_seq_orientation);
        }
        if link_type != ReadLinkType::None {
            let phase_set = self.phase_set.expect("phase set");
            for secondary_link in phase_set.secondary_links().clone() {
                if !secondary_link.has_assembly(assembly) {
                    continue;
                }
                let secondary_assembly = secondary_link.other_assembly(assembly);
                let mut link_junction_offset = 0;
                if !secondary_link.inserted_bases().is_empty() {
                    link_junction_offset += secondary_link.inserted_bases().len() as i32;
                } else if !secondary_link.overlap_bases().is_empty() {
                    link_junction_offset -= secondary_link.overlap_bases().len() as i32;
                }
                let mut secondary_assembly_junction_seq_index = assembly_junction_seq_index;
                if link_type == ReadLinkType::Left {
                    secondary_assembly_junction_seq_index += link_junction_offset;
                } else {
                    secondary_assembly_junction_seq_index -= link_junction_offset;
                }
                self.set_assembly_read_indices(secondary_assembly, is_reversed, secondary_assembly_junction_seq_index, ReadLinkType::None);
            }
        }
    }

    fn find_matching_support_read(&mut self, read: &SupP) -> Option<SupP> {
        match self.fragment_reads_map.get_mut(read.id()) {
            None => {
                self.fragment_reads_map.insert(read.id().to_string(), vec![read.clone()]);
                None
            }
            Some(fragment_reads) => {
                let matched = fragment_reads.iter().find(|x| x.first_in_pair() == read.first_in_pair()).cloned();
                fragment_reads.push(read.clone());
                matched
            }
        }
    }

    pub fn update_sequence_info(&mut self, new_sequence: String) {
        self.full_sequence_length = new_sequence.len() as i32;
        self.full_sequence = Some(new_sequence);
    }

    pub fn is_simple_indel_assembly(&self) -> bool {
        let a = self.assemblies();
        a.len() == 2 && a.iter().all(|x| x.outcome() == AssemblyOutcome::Linked) && is_assembly_indel_link(&a[0], &a[1])
    }
    pub fn is_local_indel_assembly(&self) -> bool {
        let a = self.assemblies();
        a.len() == 2 && a.iter().all(|x| x.outcome() == AssemblyOutcome::LocalIndel)
    }
    pub fn saga_match(&self) -> Option<&SagaSequenceMatch> {
        self.saga_match.as_ref()
    }
    pub fn is_saga_matched(&self) -> bool {
        self.saga_match.is_some()
    }
    pub fn set_saga_match(&mut self, m: Option<SagaSequenceMatch>) {
        self.saga_match = m;
    }
}

impl std::fmt::Display for AssemblyAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} length({}) breakends({})", self.info(), self.full_sequence_length, self.breakends.len())
    }
}

fn add_final_facing_extension_bases(assembly: AsmP, full_sequence: &mut String, sequence_cigar: &mut Vec<CigarElement>, current_seq_length: i32, assembly_reversed: bool, assembly_links: &[LinkP]) {
    if assembly_links[0].find_matched_assembly(assembly).is_some() {
        return;
    }
    let assembly_extension_bases = if assembly.is_reverse_junction() { reverse_complement_str(&assembly.form_junction_sequence()) } else { assembly.form_junction_sequence() };
    full_sequence.push_str(&assembly_extension_bases);
    log_build_info(assembly, current_seq_length, assembly_extension_bases.len() as i32, assembly_reversed, "outer-ext-bases");
    build_sequence_cigar(sequence_cigar, CigarOp::S, assembly_extension_bases.len() as i32);
}

fn build_sequence_cigar(elements: &mut Vec<CigarElement>, operator: CigarOp, length: i32) {
    if length == 0 {
        return;
    }
    match elements.last_mut() {
        Some(last) if last.op == operator => last.len += length,
        _ => elements.push(CigarElement::new(length, operator)),
    }
}

fn log_build_info(assembly: AsmP, current_seq_length: i32, assembly_base_length: i32, is_reversed: bool, other_info: &str) {
    if assembly_build_debug() {
        crate::debug!(
            "seqLength({} -> {}) adding assembly({}) {} {}",
            current_seq_length, current_seq_length + assembly_base_length, assembly.junction().coords(), if is_reversed { "rev" } else { "fwd" }, other_info
        );
    }
}
