//! PhaseSetBuilder: within a phase group, finds links between assemblies (local indels, split links, facing
//! links), extends assemblies with unmapped and remote reads, and chains the links into phase sets.

use std::cmp::Ordering;
use std::collections::HashSet;

use crate::common::compat;
use crate::assembly::assembly_utils::{self, find_matching_assembly};
use crate::assembly::common_utils::{is_discordant_fragment as record_is_discordant_fragment, is_line_insert_pair};
use crate::assembly::constants::*;
use crate::assembly::globals::{clear_cached_read, has_paired_reads, phase_processing_limit, write_candidate_reads};
use crate::assembly::phase::assembly_linker::{is_assembly_indel_link, is_facing_assembly_candidate, try_assembly_facing, try_assembly_indel, try_assembly_overlap};
use crate::assembly::phase::local_sequence_matcher::LocalSequenceMatcher;
use crate::assembly::phase::remote_read_extractor::{collect_candidate_remote_regions, purge_supplementary_reads, RemoteReadExtractor};
use crate::assembly::phase::unmapped_base_extender::UnmappedBaseExtender;
use crate::assembly::ptr::A;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::is_discordant_fragment;
use crate::assembly::ref_base_extender::{check_add_ref_base_read, extend_ref_bases, trim_assembly_ref_bases};
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::types::assembly_link::{AssemblyLink, LinkP};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::phase_group::{PhaseGroup, PhaseGroupP};
use crate::assembly::types::phase_set::{PhaseSet, PhaseSetP};
use crate::assembly::types::remote_region::RemoteRegionP;
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::{AssemblyOutcome, LinkType, StructuralVariantType, SupportType};
use crate::common::java_id_hash::IdHashSet;
use crate::common::java_sort::java_sort;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ExtensionType {
    Indel,
    LocalRefMatch,
    LocalDelDup,
    SplitLink,
    Unmapped,
}

impl ExtensionType {
    fn ordinal(self) -> i32 {
        self as i32
    }
}

pub struct ExtensionCandidate {
    pub type_: ExtensionType,
    pub assembly: AsmP,
    pub second_assembly: Option<AsmP>,
    pub link: Option<LinkP>,
    pub support_count: i32,
    pub extra_info: String,
    pub extender: Option<UnmappedBaseExtender>,
    is_valid: bool,
    selected: bool,
}

pub type ExtCandP = A<ExtensionCandidate>;

impl ExtensionCandidate {
    pub fn from_assemblies(type_: ExtensionType, assembly1: AsmP, assembly2: AsmP) -> ExtCandP {
        A::new(ExtensionCandidate { type_, assembly: assembly1, second_assembly: Some(assembly2), link: None, support_count: 0, extra_info: String::new(), extender: None, is_valid: false, selected: false })
    }
    pub fn from_link(type_: ExtensionType, link: LinkP) -> ExtCandP {
        A::new(ExtensionCandidate {
            type_,
            assembly: link.first(),
            second_assembly: Some(link.second()),
            link: Some(link),
            support_count: 0,
            extra_info: String::new(),
            extender: None,
            is_valid: true,
            selected: false,
        })
    }
    pub fn from_extender(type_: ExtensionType, assembly: AsmP, extender: UnmappedBaseExtender, support_count: i32) -> ExtCandP {
        A::new(ExtensionCandidate { type_, assembly, second_assembly: None, link: None, support_count, extra_info: String::new(), extender: Some(extender), is_valid: true, selected: false })
    }
    pub fn selected(&self) -> bool {
        self.selected
    }
    pub fn mark_selected(&mut self) {
        self.selected = true;
    }
    pub fn matches_assemblies(&self, assembly1: AsmP, assembly2: AsmP) -> bool {
        (self.assembly == assembly1 && self.second_assembly == Some(assembly2)) || (self.assembly == assembly2 && self.second_assembly == Some(assembly1))
    }
    pub fn is_valid(&self) -> bool {
        if !self.is_valid {
            return false;
        }
        self.support_count > 0 || self.type_ == ExtensionType::LocalDelDup
    }
}

/// ExtensionCandidate.StandardComparator
pub fn standard_comparator(first: &ExtCandP, second: &ExtCandP) -> Ordering {
    if first.support_count != second.support_count {
        return if first.support_count > second.support_count { Ordering::Less } else { Ordering::Greater };
    }
    // Java: the junction offset of the first candidate is compared with itself, so this falls through to the type
    if !compat::java("A19") && first.type_ == ExtensionType::SplitLink && second.type_ == ExtensionType::SplitLink {
        let junc_offset = |c: &ExtCandP| c.link.map(|l| l.inserted_bases().len().max(l.overlap_bases().len())).unwrap_or(0);
        let (first_junc_offset, second_junc_offset) = (junc_offset(first), junc_offset(second));
        if first_junc_offset != second_junc_offset {
            // favour a more precise link
            return first_junc_offset.cmp(&second_junc_offset);
        }
    }
    (-first.type_.ordinal()).cmp(&(-second.type_.ordinal()))
}

/// ExtensionCandidate.LocalLinkComparator
pub fn local_link_comparator(first: &ExtCandP, second: &ExtCandP) -> Ordering {
    if first.type_ != second.type_ {
        return first.type_.ordinal().cmp(&second.type_.ordinal());
    }
    (-first.support_count).cmp(&(-second.support_count))
}

const HIGH_ASSEMBLY_COUNT: usize = 100;
const ULTRA_HIGH_ASSEMBLY_COUNT: usize = 1000;
const HIGH_ASSEMBLY_READ_COUNT: usize = 100;
const MAX_HIGH_ASSEMBLY_MATCHED_READS: usize = 1000;

pub struct PhaseSetBuilder<'a, 'b> {
    phase_group: PhaseGroupP,
    ref_genome: Option<&'a dyn RefGenomeInterface>,
    remote_read_extractor: &'a mut RemoteReadExtractor<'b>,
    local_sequence_matcher: LocalSequenceMatcher<'a>,
    /// the top-scored assemblies of a very large group; otherwise the phase group's own (live) list is used
    filtered_assemblies: Option<Vec<AsmP>>,
    extension_candidates: Vec<ExtCandP>,
    line_related_assemblies: IdHashSet<AsmP>,
    branched_assemblies: Vec<AsmP>,
    locally_linked_assemblies: HashSet<AsmP>,
    split_links: Vec<LinkP>,
    facing_links: Vec<LinkP>,
    has_high_assembly_count: bool,
    routine_iteration: i32,
}

fn assembly_has_high_read_count(assembly: AsmP) -> bool {
    assembly.support_count() >= HIGH_ASSEMBLY_READ_COUNT || assembly.candidate_support().len() >= HIGH_ASSEMBLY_READ_COUNT
}

fn populate_read_ids(assembly: AsmP, read_ids: &mut HashSet<String>, local_indel_only: bool) {
    if assembly_has_high_read_count(assembly) {
        if local_indel_only {
            for x in assembly.support().iter().filter(|x| !x.is_discordant() && x.is_mate_mapped()) {
                read_ids.insert(x.id().to_string());
            }
            for x in assembly.candidate_support().iter().filter(|x| !is_discordant_fragment(x) && x.is_mate_mapped()) {
                read_ids.insert(x.id().to_string());
            }
        }
        // otherwise populated from the regions shared with another, non-local assembly
    } else {
        for x in assembly.support() {
            read_ids.insert(x.id().to_string());
        }
        if !local_indel_only {
            for x in assembly.candidate_support() {
                read_ids.insert(x.id().to_string());
            }
        }
    }
}

fn populate_non_local_read_ids(first_assembly: AsmP, second_assembly: AsmP, first_read_ids: &mut HashSet<String>, second_read_ids: &mut HashSet<String>) {
    let first_matched_regions: Vec<RemoteRegionP> = first_assembly.remote_regions().iter().filter(|x| x.overlaps_assembly(&second_assembly)).copied().collect();
    let second_matched_regions: Vec<RemoteRegionP> = second_assembly.remote_regions().iter().filter(|x| x.overlaps_assembly(&first_assembly)).copied().collect();
    if first_matched_regions.is_empty() || second_matched_regions.is_empty() {
        return;
    }
    for x in &first_matched_regions {
        first_read_ids.extend(x.read_ids().iter_unordered().cloned());
    }
    for x in &second_matched_regions {
        second_read_ids.extend(x.read_ids().iter_unordered().cloned());
    }
}

fn has_shared_fragment(first_read_ids: &HashSet<String>, second_read_ids: &HashSet<String>) -> bool {
    first_read_ids.iter().any(|x| second_read_ids.contains(x))
}

fn count_shared_fragments(first_read_ids: &HashSet<String>, second_read_ids: &HashSet<String>) -> i32 {
    second_read_ids.iter().filter(|x| first_read_ids.contains(*x)).count() as i32
}

pub fn is_local_assembly_candidate(first: AsmP, second: AsmP) -> bool {
    assembly_utils::is_local_assembly_candidate(&first, &second, false, false)
}

fn remove_first<T: PartialEq>(list: &mut Vec<T>, item: &T) -> bool {
    match list.iter().position(|x| x == item) {
        Some(i) => {
            list.remove(i);
            true
        }
        None => false,
    }
}

/// hmf-common SvUtils.isShortLocalDelDupIns
fn is_short_local_del_dup_ins(sv_type: StructuralVariantType, sv_length: i32) -> bool {
    matches!(sv_type, StructuralVariantType::Del | StructuralVariantType::Dup | StructuralVariantType::Ins) && sv_length <= 1000
}

impl<'a, 'b> PhaseSetBuilder<'a, 'b> {
    pub fn new(ref_genome: Option<&'a dyn RefGenomeInterface>, remote_read_extractor: &'a mut RemoteReadExtractor<'b>, phase_group: PhaseGroupP) -> Self {
        remote_read_extractor.reset_counts();
        let assembly_count = phase_group.assemblies().len();
        let has_high_assembly_count = assembly_count >= HIGH_ASSEMBLY_COUNT;
        let filtered_assemblies = if assembly_count >= ULTRA_HIGH_ASSEMBLY_COUNT { Some(filter_ultra_high_assemblies(phase_group)) } else { None };
        PhaseSetBuilder {
            phase_group,
            ref_genome,
            remote_read_extractor,
            local_sequence_matcher: LocalSequenceMatcher::new(ref_genome, LOCAL_ASSEMBLY_MATCH_DISTANCE),
            filtered_assemblies,
            extension_candidates: Vec::new(),
            line_related_assemblies: IdHashSet::new(),
            branched_assemblies: Vec::new(),
            locally_linked_assemblies: HashSet::new(),
            split_links: Vec::new(),
            facing_links: Vec::new(),
            has_high_assembly_count,
            routine_iteration: 0,
        }
    }

    /// mAssemblies: the phase group's list, which grows when assemblies are branched
    fn assemblies(&self) -> &Vec<AsmP> {
        match &self.filtered_assemblies {
            Some(a) => a,
            None => self.phase_group.assemblies(),
        }
    }
    fn assembly_count(&self) -> usize {
        self.assemblies().len()
    }
    fn assembly_at(&self, i: usize) -> AsmP {
        self.assemblies()[i]
    }

    pub fn run(&mut self) {
        self.build_phase_sets();
        PhaseGroup::finalise_phase_set_alignments(self.phase_group);
    }

    pub fn build_phase_sets(&mut self) {
        self.find_line_extensions();
        self.find_local_links();
        self.find_other_links_and_extensions();
        self.add_unlinked_assembly_ref_support();
        self.form_facing_links();
        self.check_branched_assemblies();
        self.form_phase_sets();
        self.add_chained_support();
        self.cleanup_assemblies();
    }

    fn find_local_links(&mut self) {
        self.routine_iteration = 0;
        self.find_split_link_candidates(true, None);
        java_sort(&mut self.extension_candidates, local_link_comparator);
        for extension_candidate in self.extension_candidates.clone() {
            if !extension_candidate.is_valid() {
                continue;
            }
            let assembly_link = extension_candidate.link.expect("candidate link");
            extension_candidate.m().mark_selected();
            let (a1, a2) = (extension_candidate.assembly, extension_candidate.second_assembly.expect("second assembly"));
            let is_primary_link = !self.locally_linked_assemblies.contains(&a1) && !self.locally_linked_assemblies.contains(&a2);
            self.apply_split_link(assembly_link, is_primary_link, true);
            self.locally_linked_assemblies.insert(a1);
            self.locally_linked_assemblies.insert(a2);
        }
        for assembly in self.assemblies().clone() {
            if !self.locally_linked_assemblies.contains(&assembly) {
                self.forms_local_link(assembly);
            }
        }
        self.extension_candidates.clear();
    }

    fn forms_local_link(&mut self, assembly: AsmP) -> bool {
        if assembly.discordant_only() {
            return false;
        }
        let Some(local_ref_link) = self.local_sequence_matcher.try_local_assembly_link(assembly) else { return false };
        assembly.m().set_outcome(AssemblyOutcome::LocalIndel);
        self.locally_linked_assemblies.insert(assembly);
        self.split_links.push(local_ref_link);
        let local_ref_assembly = local_ref_link.other_assembly(assembly);
        local_ref_assembly.m().set_outcome(AssemblyOutcome::LocalIndel);
        local_ref_assembly.m().set_phase_group(Some(self.phase_group));
        true
    }

    fn find_split_link_candidates(&mut self, check_local_indels: bool, extended_assemblies: Option<&HashSet<AsmP>>) {
        if self.assembly_count() < 2 {
            return;
        }
        let mut i = 0;
        while i + 1 < self.assembly_count() {
            let assembly1 = self.assembly_at(i);
            i += 1;
            if !check_local_indels && !self.allow_local_indel_secondary_links(assembly1) {
                continue;
            }
            let mut first_read_ids: HashSet<String> = HashSet::new();
            let first_high_read_count = assembly_has_high_read_count(assembly1);
            if !check_local_indels {
                populate_read_ids(assembly1, &mut first_read_ids, check_local_indels);
            }
            let mut j = i;
            while j < self.assembly_count() {
                let assembly2 = self.assembly_at(j);
                j += 1;
                if !check_local_indels && !self.allow_local_indel_secondary_links(assembly2) {
                    continue;
                }
                let mut either_extended = false;
                if let Some(extended) = extended_assemblies {
                    if !extended.contains(&assembly1) && !extended.contains(&assembly2) {
                        continue;
                    }
                    if self.split_links.iter().any(|x| x.matches(assembly1, assembly2)) {
                        continue;
                    }
                    either_extended = true;
                }
                let is_indel_junction_link = is_assembly_indel_link(&assembly1, &assembly2);
                let mut is_local_link_candidate = false;
                if !assembly1.discordant_only() && !assembly2.discordant_only() && assembly1.indel() == assembly2.indel() {
                    is_local_link_candidate = is_indel_junction_link || is_local_assembly_candidate(assembly1, assembly2);
                }
                if check_local_indels {
                    if !is_local_link_candidate {
                        continue;
                    }
                } else if is_local_link_candidate && !either_extended {
                    continue;
                }
                let mut second_read_ids: HashSet<String> = HashSet::new();
                let second_high_read_count = assembly_has_high_read_count(assembly2);
                if !check_local_indels {
                    populate_read_ids(assembly2, &mut second_read_ids, check_local_indels);
                }
                // either separate sets from overlapping remote regions, or the two read id sets themselves
                let mut related: Option<(HashSet<String>, HashSet<String>)> = None;
                if !check_local_indels && (first_high_read_count || second_high_read_count) {
                    let mut f = HashSet::new();
                    let mut s = HashSet::new();
                    populate_non_local_read_ids(assembly1, assembly2, &mut f, &mut s);
                    related = Some((f, s));
                }
                let requires_shared_fragments = !check_local_indels && !(either_extended && is_local_link_candidate);
                if requires_shared_fragments {
                    let shared = match &related {
                        Some((f, s)) => has_shared_fragment(f, s),
                        None => has_shared_fragment(&first_read_ids, &second_read_ids),
                    };
                    if !shared {
                        continue;
                    }
                }
                let Some(assembly_link) = self.check_split_link(assembly1, assembly2, is_local_link_candidate) else {
                    if check_local_indels {
                        self.extension_candidates.push(ExtensionCandidate::from_assemblies(ExtensionType::LocalDelDup, assembly1, assembly2));
                    }
                    continue;
                };
                let type_ = if is_local_link_candidate {
                    if is_indel_junction_link { ExtensionType::Indel } else { ExtensionType::LocalDelDup }
                } else {
                    ExtensionType::SplitLink
                };
                let extension_candidate = ExtensionCandidate::from_link(type_, assembly_link);
                self.extension_candidates.push(extension_candidate.clone());
                if check_local_indels {
                    if first_read_ids.is_empty() {
                        populate_read_ids(assembly1, &mut first_read_ids, check_local_indels);
                    }
                    populate_read_ids(assembly2, &mut second_read_ids, check_local_indels);
                }
                let count = match &related {
                    Some((f, s)) => count_shared_fragments(f, s),
                    None => count_shared_fragments(&first_read_ids, &second_read_ids),
                };
                extension_candidate.m().support_count += count;
            }
            if self.check_max_routine_iteration() {
                return;
            }
        }
    }

    fn allow_local_indel_secondary_links(&self, assembly: AsmP) -> bool {
        if self.line_related_assemblies.contains(&assembly) {
            return true;
        }
        if !self.locally_linked_assemblies.contains(&assembly) {
            return true;
        }
        match self.split_links.iter().find(|x| x.has_assembly(assembly)) {
            None => true,
            Some(existing_link) => existing_link.inserted_bases().len() as i32 > existing_link.length(),
        }
    }

    fn find_other_links_and_extensions(&mut self) {
        self.find_unmapped_extensions();
        let extended_assemblies = self.apply_other_links_and_extensions(None);
        if !extended_assemblies.is_empty() {
            self.apply_other_links_and_extensions(Some(&extended_assemblies));
        }
    }

    fn apply_other_links_and_extensions(&mut self, extended_assemblies: Option<&HashSet<AsmP>>) -> HashSet<AsmP> {
        self.routine_iteration = 0;
        self.find_split_link_candidates(false, extended_assemblies);
        self.routine_iteration = 0;
        let check_only_extended_assemblies = extended_assemblies.is_some();
        let mut remaining_candidates: Vec<ExtCandP> = self
            .extension_candidates
            .iter()
            .filter(|x| !x.selected())
            .filter(|x| x.is_valid())
            .filter(|x| !check_only_extended_assemblies || x.type_ != ExtensionType::Unmapped)
            .cloned()
            .collect();
        if remaining_candidates.is_empty() {
            return HashSet::new();
        }
        java_sort(&mut remaining_candidates, standard_comparator);
        let mut primary_linked_assemblies: HashSet<AsmP> = self.locally_linked_assemblies.clone();
        for x in &self.split_links {
            primary_linked_assemblies.insert(x.first());
            primary_linked_assemblies.insert(x.second());
        }
        let mut newly_extended_assemblies: HashSet<AsmP> = HashSet::new();
        for extension_candidate in remaining_candidates {
            match extension_candidate.type_ {
                ExtensionType::SplitLink => {
                    let (a1, a2) = (extension_candidate.assembly, extension_candidate.second_assembly.expect("second assembly"));
                    if let Some(extended) = extended_assemblies {
                        if !extended.contains(&a1) && !extended.contains(&a2) {
                            continue;
                        }
                    }
                    let either_in_primary = primary_linked_assemblies.contains(&a1) || primary_linked_assemblies.contains(&a2);
                    if either_in_primary && (a1.discordant_only() || a2.discordant_only()) {
                        continue;
                    }
                    extension_candidate.m().mark_selected();
                    self.apply_split_link(extension_candidate.link.expect("candidate link"), !either_in_primary, false);
                    if !either_in_primary {
                        primary_linked_assemblies.insert(a1);
                        primary_linked_assemblies.insert(a2);
                    }
                }
                ExtensionType::Unmapped => {
                    if extension_candidate.assembly.outcome() == AssemblyOutcome::Unset {
                        extension_candidate.m().mark_selected();
                        self.apply_unmapped_read_extension(&extension_candidate);
                        newly_extended_assemblies.insert(extension_candidate.assembly);
                    }
                }
                ExtensionType::LocalDelDup => {
                    extension_candidate.m().mark_selected();
                    let (a1, a2) = (extension_candidate.assembly, extension_candidate.second_assembly.expect("second assembly"));
                    let is_primary_link = !self.locally_linked_assemblies.contains(&a1) && !self.locally_linked_assemblies.contains(&a2);
                    self.apply_split_link(extension_candidate.link.expect("candidate link"), is_primary_link, true);
                    self.locally_linked_assemblies.insert(a1);
                    self.locally_linked_assemblies.insert(a2);
                }
                _ => {}
            }
        }
        newly_extended_assemblies
    }

    fn find_unmapped_extensions(&mut self) {
        if !has_paired_reads() {
            return;
        }
        self.routine_iteration = 0;
        let java_a20 = compat::java("A20");
        for assembly in self.assemblies().clone() {
            if self.locally_linked_assemblies.contains(&assembly) || self.line_related_assemblies.contains(&assembly) {
                continue;
            }
            let mut unmapped_reads: Vec<ReadP> = assembly.unmapped_reads().clone();
            let mut remote_regions = collect_candidate_remote_regions(assembly, self.assemblies(), self.has_high_assembly_count);
            let mut remote_reads = self.remote_read_extractor.extract_remote_region_reads(self.phase_group.id(), &mut remote_regions, self.has_high_assembly_count);
            purge_supplementary_reads(&assembly, &mut remote_reads);
            unmapped_reads.extend(remote_reads);
            if unmapped_reads.is_empty() {
                // Java: the first assembly without any reads ends the search for all the assemblies after it
                if java_a20 {
                    break;
                }
                continue;
            }
            let mut unmapped_base_extender = UnmappedBaseExtender::new(assembly);
            unmapped_base_extender.process_reads(&mut unmapped_reads);
            if !unmapped_base_extender.support_reads().is_empty() {
                let support = unmapped_base_extender.support_reads().len() as i32;
                let read_span = unmapped_base_extender.extension_bases().len();
                let extension_candidate = ExtensionCandidate::from_extender(ExtensionType::Unmapped, assembly, unmapped_base_extender, support);
                extension_candidate.m().extra_info = format!("readSpan({read_span})");
                self.extension_candidates.push(extension_candidate);
            }
            if self.check_max_routine_iteration() {
                break;
            }
        }
    }

    fn find_line_extensions(&mut self) {
        for x in self.assemblies().clone() {
            if x.has_line_sequence() {
                self.line_related_assemblies.insert(x);
            }
        }
        if self.line_related_assemblies.is_empty() {
            return;
        }
        self.routine_iteration = 0;
        let mut proximate_line_assemblies: Vec<AsmP> = Vec::new();
        let mut extended_assemblies: HashSet<AsmP> = HashSet::new();
        for assembly in self.assemblies().clone() {
            if assembly.has_line_sequence() {
                continue;
            }
            let Some(line_assembly) = self.line_related_assemblies.iter().find(|x| is_line_insert_pair(&assembly, x)).copied() else { continue };
            proximate_line_assemblies.push(assembly);
            if !has_paired_reads() {
                continue;
            }
            let mut shared_unmapped_reads: Vec<ReadP> = Vec::new();
            let mut combined_remote_regions: Vec<RemoteRegionP> = Vec::new();
            let line_site_assemblies = [assembly, line_assembly];
            for line_site_assembly in line_site_assemblies {
                shared_unmapped_reads.extend(line_site_assembly.unmapped_reads().iter().cloned());
                combined_remote_regions.extend(line_site_assembly.remote_regions().iter().filter(|x| !x.is_supp_only_region()).copied());
            }
            let mut remote_reads = self.remote_read_extractor.extract_remote_region_reads(self.phase_group.id(), &mut combined_remote_regions, self.has_high_assembly_count);
            purge_supplementary_reads(&assembly, &mut remote_reads);
            shared_unmapped_reads.extend(remote_reads);
            if shared_unmapped_reads.is_empty() {
                continue;
            }
            for line_site_assembly in line_site_assemblies {
                if !extended_assemblies.insert(line_site_assembly) {
                    continue;
                }
                let mut unmapped_base_extender = UnmappedBaseExtender::new(line_site_assembly);
                let mut reads_copy = shared_unmapped_reads.clone();
                unmapped_base_extender.process_reads(&mut reads_copy);
                if !unmapped_base_extender.support_reads().is_empty() {
                    line_site_assembly.m().expand_extension_bases(&unmapped_base_extender.extension_bases(), &unmapped_base_extender.base_qualities(), unmapped_base_extender.support_reads());
                }
            }
        }
        for a in proximate_line_assemblies {
            self.line_related_assemblies.insert(a);
        }
    }

    fn apply_split_link(&mut self, assembly_link: LinkP, is_primary_link: bool, is_local_indel: bool) {
        let allow_branching = is_primary_link && !(assembly_link.sv_type() == StructuralVariantType::Dup && assembly_link.length() < PROXIMATE_DUP_LENGTH);
        let allow_discordant_reads = !is_local_indel && !is_short_local_del_dup_ins(assembly_link.sv_type(), assembly_link.length());
        if is_local_indel {
            self.adjust_local_indel_support(assembly_link.first(), assembly_link.second());
        }
        self.apply_split_link_support(assembly_link.first(), assembly_link.second(), allow_branching, allow_discordant_reads);
        if is_primary_link {
            self.split_links.push(assembly_link);
            assembly_link.first().m().set_outcome(AssemblyOutcome::Linked);
            assembly_link.second().m().set_outcome(AssemblyOutcome::Linked);
        } else {
            self.phase_group.m().secondary_links_mut().push(assembly_link);
            assembly_link.first().m().set_outcome(AssemblyOutcome::Secondary);
            assembly_link.second().m().set_outcome(AssemblyOutcome::Secondary);
        }
    }

    fn apply_unmapped_read_extension(&mut self, extension_candidate: &ExtCandP) {
        let unmapped_base_extender = extension_candidate.extender.as_ref().expect("unmapped base extender");
        let assembly = extension_candidate.assembly;
        if !unmapped_base_extender.support_reads().is_empty() {
            assembly.m().expand_extension_bases(&unmapped_base_extender.extension_bases(), &unmapped_base_extender.base_qualities(), unmapped_base_extender.support_reads());
        }
    }

    fn check_split_link(&self, assembly1: AsmP, assembly2: AsmP, is_local_indel: bool) -> Option<LinkP> {
        if assembly1.junction() == assembly2.junction() {
            return None;
        }
        if is_local_indel {
            if let Some(link) = try_assembly_indel(assembly1, assembly2) {
                return Some(link);
            }
        }
        try_assembly_overlap(assembly1, assembly2, true, is_local_indel)
    }

    fn add_unlinked_assembly_ref_support(&mut self) {
        self.routine_iteration = 0;
        let assemblies = self.assemblies().clone();
        let allow_ref_side_soft_clip_branching = true;
        for assembly in assemblies {
            if assembly.outcome() != AssemblyOutcome::Unset {
                continue;
            }
            if self.split_links.iter().any(|x| x.has_assembly(assembly)) {
                continue;
            }
            let mut ref_extension_reads: Vec<ReadP> = Vec::new();
            let extension_support: Vec<SupP> = assembly.support().iter().filter(|x| x.type_() == SupportType::Extension).cloned().collect();
            for read in assembly.candidate_support() {
                if read.has_junction_mate() || extension_support.iter().any(|x| x.matches_fragment_read(read, false)) {
                    ref_extension_reads.push(read.clone());
                }
            }
            extend_ref_bases(assembly, &mut ref_extension_reads, self.ref_genome, allow_ref_side_soft_clip_branching);
        }
    }

    fn apply_split_link_support(&mut self, assembly1: AsmP, assembly2: AsmP, allow_branching: bool, allow_discordant_reads: bool) -> bool {
        if !has_paired_reads() {
            return false;
        }
        let mut matched_candidates1: Vec<ReadP> = Vec::new();
        let mut matched_candidates2: Vec<ReadP> = Vec::new();
        let has_high_assembly_count = self.has_high_assembly_count;
        let candidate_support_limited1 = has_high_assembly_count && assembly1.candidate_support().len() > MAX_HIGH_ASSEMBLY_MATCHED_READS;
        let candidate_support_limited2 = has_high_assembly_count && assembly2.candidate_support().len() > MAX_HIGH_ASSEMBLY_MATCHED_READS;
        // the Java sub-lists are windows over the first reads of the assemblies' own lists
        let mut window1 = if candidate_support_limited1 { MAX_HIGH_ASSEMBLY_MATCHED_READS } else { assembly1.candidate_support().len() };
        let mut window2 = if candidate_support_limited2 { MAX_HIGH_ASSEMBLY_MATCHED_READS } else { assembly2.candidate_support().len() };
        check_matching_candidate_support(assembly2, allow_discordant_reads, assembly1, &mut window1, Some((assembly2, &mut window2)), &mut matched_candidates1, &mut matched_candidates2);
        check_matching_candidate_support(assembly1, allow_discordant_reads, assembly2, &mut window2, None, &mut matched_candidates2, &mut matched_candidates1);
        // Java: takes the matched reads out of the assemblies' lists a second time, on the assumption that the
        // sub-lists were copies; they are views, so the reads are gone already and this removes a second
        // occurrence of the same read, if the list holds one
        if compat::java("A21") {
            if candidate_support_limited1 {
                for x in &matched_candidates1 {
                    remove_first(assembly1.m().candidate_support_mut(), x);
                }
            }
            if candidate_support_limited2 {
                for x in &matched_candidates2 {
                    remove_first(assembly2.m().candidate_support_mut(), x);
                }
            }
        }
        add_matching_extension_candidates(assembly1, &mut matched_candidates1);
        add_matching_extension_candidates(assembly2, &mut matched_candidates2);

        let mut ref_candidates1: Vec<ReadP> = Vec::new();
        let mut has_non_local_tumor_fragment = false;
        let mut has_non_local_ref_fragment = false;
        for read in &matched_candidates1 {
            let disc = record_is_discordant_fragment(read.bam_record(), -1, read.supplementary_data());
            if read.is_reference() {
                ref_candidates1.push(read.clone());
                has_non_local_ref_fragment |= disc;
            } else {
                has_non_local_tumor_fragment |= disc;
            }
        }
        if has_non_local_tumor_fragment && !has_non_local_ref_fragment {
            let ref_candidates2: Vec<ReadP> = matched_candidates2.iter().filter(|x| x.is_reference()).cloned().collect();
            for x in &ref_candidates1 {
                remove_first(&mut matched_candidates1, x);
            }
            for x in &ref_candidates2 {
                remove_first(&mut matched_candidates2, x);
            }
        }
        if matched_candidates1.is_empty() && matched_candidates2.is_empty() {
            return false;
        }
        extend_ref_bases(assembly1, &mut matched_candidates1, self.ref_genome, allow_branching);
        extend_ref_bases(assembly2, &mut matched_candidates2, self.ref_genome, allow_branching);
        for assembly in self.assemblies().clone() {
            if assembly.outcome() == AssemblyOutcome::DupBranched && !self.branched_assemblies.contains(&assembly) {
                self.branched_assemblies.push(assembly);
            }
        }
        true
    }

    fn adjust_local_indel_support(&self, assembly1: AsmP, assembly2: AsmP) {
        if !is_local_assembly_candidate(assembly1, assembly2) {
            return;
        }
        for is_start in [true, false] {
            let assembly = if is_start { assembly1 } else { assembly2 };
            let other_assembly = if is_start { assembly2 } else { assembly1 };
            for mate_read in assembly.concordant_candidates().clone() {
                if other_assembly.is_forward_junction() {
                    if mate_read.orientation().is_forward() && mate_read.alignment_end() <= other_assembly.junction().position {
                        other_assembly.m().add_candidate_support(mate_read);
                    }
                } else if mate_read.orientation().is_reverse() && mate_read.alignment_start() >= other_assembly.junction().position {
                    other_assembly.m().add_candidate_support(mate_read);
                }
            }
        }
    }

    fn form_facing_links(&mut self) {
        if self.assembly_count() == 1 || (self.assembly_count() == 2 && self.split_links.len() == 1) {
            return;
        }
        let mut facing_assemblies: HashSet<AsmP> = HashSet::new();
        let mut split_links: Vec<LinkP> = self.split_links.clone();
        split_links.extend(self.phase_group.secondary_links().iter().copied());
        let mut i = 0;
        while i + 1 < self.assembly_count() {
            let assembly1 = self.assembly_at(i);
            i += 1;
            if !is_facing_assembly_candidate(assembly1, facing_assemblies.contains(&assembly1), &self.split_links) {
                continue;
            }
            let mut j = i;
            while j < self.assembly_count() {
                let assembly2 = self.assembly_at(j);
                j += 1;
                if !is_facing_assembly_candidate(assembly2, facing_assemblies.contains(&assembly2), &self.split_links) {
                    continue;
                }
                let Some(facing_link) = try_assembly_facing(assembly1, assembly2, &split_links) else { continue };
                self.facing_links.push(facing_link);
                facing_assemblies.insert(assembly1);
                facing_assemblies.insert(assembly2);
            }
        }
    }

    fn check_branched_assemblies(&mut self) {
        for assembly in self.branched_assemblies.clone() {
            let has_split_link = self.split_links.iter().any(|x| x.has_assembly(assembly));
            let has_facing_link = self.facing_links.iter().any(|x| x.has_assembly(assembly));
            if !has_split_link && !has_facing_link {
                continue;
            }
            if has_facing_link && !has_split_link {
                let Some(original_assembly) = find_matching_assembly(self.assemblies(), assembly, true) else { continue };
                let existing_link = self.split_links.iter().find(|x| x.has_assembly(original_assembly)).copied();
                let orig_has_facing_link = self.facing_links.iter().any(|x| x.has_assembly(original_assembly));
                if let (false, Some(existing_link)) = (orig_has_facing_link, existing_link) {
                    let replacement_link = AssemblyLink::swap_assemblies(existing_link, original_assembly, assembly);
                    remove_first(&mut self.split_links, &existing_link);
                    self.split_links.push(replacement_link);
                    remove_first(self.phase_group.m().derived_assemblies_mut(), &assembly);
                    original_assembly.m().set_outcome_override(AssemblyOutcome::DupBranched, true);
                    assembly.m().set_outcome(AssemblyOutcome::Linked);
                }
            }
        }
    }

    fn form_phase_sets(&mut self) {
        while !self.split_links.is_empty() || !self.facing_links.is_empty() {
            let assembly_link = if !self.split_links.is_empty() { self.split_links.remove(0) } else { self.facing_links.remove(0) };
            let phase_set = PhaseSet::new(assembly_link);
            self.phase_group.m().phase_sets_mut().push(phase_set);
            if assembly_link.type_() == LinkType::Split && assembly_link.first().outcome() == AssemblyOutcome::LocalIndel {
                continue;
            }
            for is_start in [true, false] {
                let mut linking_assembly = if is_start { assembly_link.first() } else { assembly_link.second() };
                let mut find_split = assembly_link.type_() == LinkType::Facing;
                loop {
                    let Some(next_link) = self.find_linked_assembly(linking_assembly, find_split) else { break };
                    if is_start {
                        phase_set.m().add_assembly_link_start(next_link);
                    } else {
                        phase_set.m().add_assembly_link_end(next_link);
                    }
                    find_split = !find_split;
                    linking_assembly = next_link.other_assembly(linking_assembly);
                }
            }
        }
        let mut processed_secondaries: HashSet<AsmP> = HashSet::new();
        let mut secondary_line_phase_sets: Vec<PhaseSetP> = Vec::new();
        for link in self.phase_group.secondary_links().clone() {
            let first_is_line_site = link.first().outcome() == AssemblyOutcome::Secondary && self.line_related_assemblies.contains(&link.first());
            let second_is_line_site = link.second().outcome() == AssemblyOutcome::Secondary && self.line_related_assemblies.contains(&link.second());
            if first_is_line_site || second_is_line_site {
                let processed = (first_is_line_site && processed_secondaries.contains(&link.first())) || (second_is_line_site && processed_secondaries.contains(&link.second()));
                if !processed {
                    secondary_line_phase_sets.push(PhaseSet::new(link));
                    if first_is_line_site {
                        processed_secondaries.insert(link.first());
                    }
                    if second_is_line_site {
                        processed_secondaries.insert(link.second());
                    }
                }
            } else {
                for phase_set in self.phase_group.phase_sets().clone() {
                    if phase_set.has_assembly(link.first()) || phase_set.has_assembly(link.second()) {
                        phase_set.m().add_secondary_link(link);
                    }
                }
            }
        }
        self.phase_group.m().phase_sets_mut().extend(secondary_line_phase_sets);
    }

    fn find_linked_assembly(&mut self, assembly: AsmP, find_split: bool) -> Option<LinkP> {
        let search_links = if find_split { &mut self.split_links } else { &mut self.facing_links };
        let mut index = 0;
        while index < search_links.len() {
            let link = search_links[index];
            if link.has_assembly(assembly) {
                search_links.remove(index);
                if !find_split {
                    let other_assembly = link.other_assembly(assembly);
                    let mut other_index = 0;
                    while other_index < search_links.len() {
                        let other_link = search_links[other_index];
                        if other_link.has_assembly(assembly) || other_link.has_assembly(other_assembly) {
                            search_links.remove(other_index);
                        } else {
                            other_index += 1;
                        }
                    }
                }
                return Some(link);
            }
            index += 1;
        }
        None
    }

    fn add_chained_support(&mut self) {
        for phase_set in self.phase_group.phase_sets().clone() {
            let n = phase_set.assemblies().len();
            if n <= 2 {
                continue;
            }
            for i in 0..n - 1 {
                let assembly1 = phase_set.assemblies()[i];
                let assembly_links = phase_set.find_assembly_links(assembly1);
                for j in i + 1..n {
                    let assembly2 = phase_set.assemblies()[j];
                    if assembly_links.iter().any(|x| x.has_assembly(assembly2)) {
                        continue;
                    }
                    add_matching_candidate_support(phase_set, assembly1, assembly2);
                }
            }
        }
    }

    fn cleanup_assemblies(&mut self) {
        let mut branched_assemblies_to_remove: Vec<AsmP> = Vec::new();
        for assembly in self.assemblies().clone() {
            if !write_candidate_reads() {
                assembly.m().clear_candidate_support();
            }
            assembly.m().clear_support_cached_reads();
            let mut in_phase_set = false;
            let mut in_facing_link = false;
            for phase_set in self.phase_group.phase_sets() {
                if phase_set.has_assembly(assembly) {
                    in_phase_set = true;
                    in_facing_link |= phase_set.assembly_links().iter().any(|x| x.type_() == LinkType::Facing && x.has_assembly(assembly));
                }
            }
            if in_phase_set && assembly.outcome() == AssemblyOutcome::Unset {
                assembly.m().set_outcome(AssemblyOutcome::Linked);
            }
            if !in_facing_link {
                trim_assembly_ref_bases(assembly, ASSEMBLY_REF_BASE_MAX_GAP);
            }
            if assembly.outcome() == AssemblyOutcome::DupBranched {
                if in_phase_set && in_facing_link {
                    if let Some(original_assembly) = find_matching_assembly(self.assemblies(), assembly, false) {
                        assembly.m().set_outcome(original_assembly.outcome());
                    }
                }
                if !in_facing_link {
                    branched_assemblies_to_remove.push(assembly);
                }
            }
        }
        // not in Java: with the candidate and cached reads of these assemblies dropped, their junction groups
        // can free the reads (done after the loop so that no assembly of this phase group still holds one)
        if !write_candidate_reads() && clear_cached_read() {
            for assembly in self.assemblies().clone() {
                assembly.m().candidate_reads_dropped();
            }
        }
        for branched_assembly in branched_assemblies_to_remove {
            remove_first(self.phase_group.m().assemblies_mut(), &branched_assembly);
            remove_first(self.phase_group.m().derived_assemblies_mut(), &branched_assembly);
        }
    }

    fn check_max_routine_iteration(&mut self) -> bool {
        let limit = phase_processing_limit();
        if limit == 0 {
            return false;
        }
        self.routine_iteration += 1;
        self.routine_iteration >= limit
    }
}

/// PhaseSetBuilder.checkMatchingCandidateSupport. `window` is the number of leading candidate reads of
/// `candidate_assembly` that are considered (the Java sub-list), reduced as reads are taken out of it.
fn check_matching_candidate_support(
    other_assembly: AsmP, allow_discordant_reads: bool, candidate_assembly: AsmP, window: &mut usize, mut other_candidates: Option<(AsmP, &mut usize)>,
    matched_candidates: &mut Vec<ReadP>, other_matched_candidates: &mut Vec<ReadP>,
) {
    let mut index = 0;
    while index < *window {
        let candidate_read = candidate_assembly.candidate_support()[index].clone();
        if candidate_read.has_junction_mate() {
            candidate_assembly.m().candidate_support_mut().remove(index);
            *window -= 1;
            matched_candidates.push(candidate_read);
            continue;
        }
        if !allow_discordant_reads {
            index += 1;
            continue;
        }
        if SupportRead::has_fragment_other_read(other_assembly.support(), &candidate_read) {
            candidate_assembly.m().candidate_support_mut().remove(index);
            *window -= 1;
            matched_candidates.push(candidate_read);
            continue;
        }
        if let Some((other_candidate_assembly, other_window)) = other_candidates.as_mut() {
            if **other_window > 0 {
                let matched_candidate_support: Vec<ReadP> = Read::find_matching_fragment_support(&other_candidate_assembly.candidate_support()[..**other_window], &candidate_read);
                if !matched_candidate_support.is_empty() {
                    candidate_assembly.m().candidate_support_mut().remove(index);
                    *window -= 1;
                    matched_candidates.push(candidate_read);
                    for x in &matched_candidate_support {
                        let list = other_candidate_assembly.m().candidate_support_mut();
                        if let Some(pos) = list[..**other_window].iter().position(|r| r == x) {
                            list.remove(pos);
                            **other_window -= 1;
                        }
                    }
                    other_matched_candidates.extend(matched_candidate_support);
                    continue;
                }
            }
        }
        index += 1;
    }
}

fn add_matching_extension_candidates(assembly: AsmP, matched_candidates: &mut Vec<ReadP>) {
    let mut extension_support: HashSet<String> = assembly.support().iter().filter(|x| x.type_() == SupportType::Extension).map(|x| x.id().to_string()).collect();
    if extension_support.is_empty() {
        return;
    }
    let mut index = 0;
    while index < assembly.candidate_support().len() {
        let candidate_read = assembly.candidate_support()[index].clone();
        if extension_support.remove(candidate_read.id()) {
            matched_candidates.push(candidate_read);
            assembly.m().candidate_support_mut().remove(index);
            continue;
        }
        index += 1;
    }
}

fn add_matching_candidate_support(phase_set: PhaseSetP, assembly1: AsmP, assembly2: AsmP) {
    if !phase_set.assemblies_face_in_phase_set(assembly1, assembly2) {
        return;
    }
    for is_start in [true, false] {
        let assembly = if is_start { assembly1 } else { assembly2 };
        let other_assembly = if is_start { assembly2 } else { assembly1 };
        let mut index = 0;
        while index < assembly.candidate_support().len() {
            let candidate_read = assembly.candidate_support()[index].clone();
            if candidate_read.has_junction_mate() {
                index += 1;
                continue;
            }
            let matched_read = other_assembly.support().iter().any(|x| x.matches_fragment_read(&candidate_read, false));
            if matched_read {
                assembly.m().candidate_support_mut().remove(index);
                check_add_ref_base_read(assembly, &candidate_read, SupportType::Discordant);
                continue;
            }
            if is_start {
                let matched_candidate = other_assembly.candidate_support().iter().find(|x| x.matches_fragment(&candidate_read, false)).cloned();
                if let Some(matched_candidate) = matched_candidate {
                    assembly.m().candidate_support_mut().remove(index);
                    check_add_ref_base_read(assembly, &candidate_read, SupportType::Discordant);
                    remove_first(other_assembly.m().candidate_support_mut(), &matched_candidate);
                    check_add_ref_base_read(other_assembly, &matched_candidate, SupportType::Discordant);
                    continue;
                }
            }
            index += 1;
        }
    }
}

fn score_junction_assembly(assembly: AsmP) -> f64 {
    let frag_count = assembly.support_count() as f64;
    let mismatch_read_count = assembly.mismatch_read_count();
    let ext_base_length = assembly.extension_length();
    let ref_base_length = assembly.ref_base_length();
    let frag_score = frag_count / mismatch_read_count.max(1) as f64;
    let length_score = (ext_base_length * ref_base_length) as f64 / (MIN_VARIANT_LENGTH * MIN_VARIANT_LENGTH) as f64;
    frag_score * length_score
}

/// Integers.median
fn median(values: &mut Vec<i32>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort();
    let index = values.len() / 2;
    if values.len() % 2 == 0 { (values[index - 1] + values[index]) as f64 * 0.5 } else { values[index] as f64 }
}

fn filter_ultra_high_assemblies(phase_group: PhaseGroupP) -> Vec<AsmP> {
    let mut line_count = 0;
    let mut split_count = 0;
    let mut indel_count = 0;
    let mut discordant_count = 0;
    let mut extension_lengths: Vec<i32> = Vec::new();
    let mut fragment_counts: Vec<i32> = Vec::new();
    let mut mismatch_read_counts: Vec<i32> = Vec::new();
    for assembly in phase_group.assemblies() {
        if assembly.junction().discordant_only {
            discordant_count += 1;
        } else if assembly.junction().indel_based() {
            indel_count += 1;
        } else {
            split_count += 1;
            if assembly.has_line_sequence() {
                line_count += 1;
            }
        }
        fragment_counts.push(assembly.support_count() as i32);
        extension_lengths.push(assembly.extension_length());
        mismatch_read_counts.push(assembly.mismatch_read_count());
    }
    let round = crate::assembly::types::assembly_stats::java_round;
    crate::info!(
        "pgId({}) assemblies({}) types(disc={} indel={} split={} line={}) median(frags={} extLength={} mismatchReads={})",
        phase_group.id(), phase_group.assemblies().len(), discordant_count, indel_count, split_count, line_count,
        round(median(&mut fragment_counts)) as i32, round(median(&mut extension_lengths)) as i32, round(median(&mut mismatch_read_counts)) as i32
    );
    let mut all_assemblies: Vec<AsmP> = phase_group.assemblies().clone();
    java_sort(&mut all_assemblies, |a, b| (-score_junction_assembly(*a)).partial_cmp(&(-score_junction_assembly(*b))).unwrap_or(Ordering::Equal));
    all_assemblies.truncate(ULTRA_HIGH_ASSEMBLY_COUNT);
    all_assemblies
}
