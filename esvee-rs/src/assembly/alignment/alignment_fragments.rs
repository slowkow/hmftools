//! AlignmentFragments: assigns the reads of an assembly alignment's assemblies to its breakends as split or
//! discordant support, and infers fragment lengths.

use crate::common::compat;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend::{Breakend, BreakendP, BreakendSupport};
use crate::assembly::common_utils::{position_within, positions_within};
use crate::assembly::constants::*;
use crate::assembly::globals::max_observed_concordant_frag_length;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::{AssemblyOutcome, StructuralVariantType, SupportType};
use crate::common::java_hash::JavaHashMap;
use crate::common::java_id_hash::IdHashSet;

const INVALID_DISCORANT_DISTANCE: i32 = -1;

pub struct AlignmentFragments {
    assembly_alignment: AsmAlignP,
}

struct FragmentReads {
    primary_reads: Vec<SupP>,
    duplicate_reads: Vec<SupP>,
    fragment_length: i32,
    processed: bool,
    in_primary_assembly: bool,
}

impl FragmentReads {
    fn new(read: SupP, in_primary_assembly: bool) -> Self {
        FragmentReads { primary_reads: vec![read], duplicate_reads: Vec::new(), fragment_length: -1, processed: false, in_primary_assembly }
    }
    fn add_read(&mut self, read: SupP) {
        if self.primary_reads.iter().any(|x| x.flags() == read.flags()) {
            self.duplicate_reads.push(read);
        } else {
            self.primary_reads.push(read);
        }
    }
    fn set_read_details(&self, read: &SupP) {
        read.m().set_inferred_fragment_length(self.fragment_length);
        if let Some(matching_read) = self.primary_reads.iter().find(|x| x.flags() == read.flags()) {
            let t = matching_read.breakend_support_type();
            read.m().set_breakend_support_type(t);
        }
    }
    fn set_fragment_length(&mut self) {
        self.fragment_length = self.primary_reads.iter().map(|x| x.inferred_fragment_length()).max().unwrap_or(-1);
        for x in self.primary_reads.iter().chain(self.duplicate_reads.iter()) {
            x.m().set_inferred_fragment_length(self.fragment_length);
        }
    }
    fn find_specific_read(&self, get_first: bool, use_supplementaries: bool) -> Option<SupP> {
        if let Some(read) = self.primary_reads.iter().find(|x| x.first_in_pair() == get_first && !x.is_supplementary()) {
            return Some(read.clone());
        }
        if !use_supplementaries {
            return None;
        }
        self.primary_reads.iter().find(|x| x.first_in_pair() == get_first).cloned()
    }
}

#[derive(Clone)]
struct ReadBreakendMatch {
    is_split: bool,
    breakend: BreakendP,
}

pub fn read_alignment_spans_junction(breakend: &Breakend, read: &SupportRead) -> bool {
    breakend.chromosome == read.chromosome() && read.unclipped_start() < breakend.position && read.unclipped_end() > breakend.position
}

fn is_indel(sv_type: StructuralVariantType) -> bool {
    matches!(sv_type, StructuralVariantType::Del | StructuralVariantType::Dup | StructuralVariantType::Ins)
}

fn can_use_other_read_breakend_matches(first_read: &SupP, second_read: &SupP) -> bool {
    if !first_read.is_discordant() && first_read.type_() == SupportType::Junction && second_read.type_() == SupportType::JunctionMate {
        return true;
    }
    if !second_read.is_discordant() && second_read.type_() == SupportType::Junction && first_read.type_() == SupportType::JunctionMate {
        return true;
    }
    first_read.type_() == SupportType::Extension || second_read.type_() == SupportType::Extension
}

fn add_unique_breakends(breakends: &mut IdHashSet<BreakendP>, read_breakend_matches: &[ReadBreakendMatch]) {
    for breakend_match in read_breakend_matches {
        let b = breakend_match.breakend;
        if breakends.contains(&b) {
            continue;
        }
        if let Some(other) = b.other_breakend() {
            if breakends.contains(&other) {
                continue;
            }
        }
        breakends.insert(b);
    }
}

fn read_discordant_breakend_distance(breakend: &Breakend, read: &SupportRead) -> i32 {
    if breakend.chromosome == read.chromosome() && read.orientation() == breakend.orient {
        if breakend.orient.is_forward() {
            let max_position = breakend.max_position().max(breakend.position);
            if read.alignment_end() <= max_position && read.alignment_start() >= breakend.position - DEFAULT_MAX_CONCORDANT_FRAG_LENGTH {
                return max_position - read.alignment_end();
            }
        } else {
            let min_position = breakend.min_position().min(breakend.position);
            if read.alignment_start() >= min_position && read.alignment_end() <= breakend.position + DEFAULT_MAX_CONCORDANT_FRAG_LENGTH {
                return read.alignment_start() - min_position;
            }
        }
    }
    let read_seq_index_start = read.full_assembly_index_start();
    let read_seq_index_end = read.full_assembly_index_end();
    for segment in breakend.segments() {
        if read.type_() != SupportType::Extension && !positions_within(read_seq_index_start, read_seq_index_end, segment.alignment.sequence_start(), segment.alignment.sequence_end()) {
            continue;
        }
        let mut discordant_distance = if read.full_assembly_orientation().is_forward() { segment.alignment.sequence_end() - read_seq_index_end } else { read_seq_index_start - segment.alignment.sequence_start() };
        if read.type_() == SupportType::Extension {
            discordant_distance = discordant_distance.abs();
        }
        if discordant_distance > 0 && discordant_distance < DEFAULT_MAX_CONCORDANT_FRAG_LENGTH {
            return discordant_distance;
        }
    }
    INVALID_DISCORANT_DISTANCE
}

impl AlignmentFragments {
    pub fn new(assembly_alignment: AsmAlignP, combined_sample_ids: &[String]) -> AlignmentFragments {
        for breakend in assembly_alignment.breakends() {
            let support = breakend.m().sample_support_mut();
            for _ in combined_sample_ids {
                support.push(BreakendSupport::default());
            }
        }
        AlignmentFragments { assembly_alignment }
    }

    pub fn allocate_breakend_support(&self) {
        let aa = self.assembly_alignment;
        let mut fragments: Vec<FragmentReads> = Vec::new();
        let mut fragment_map: JavaHashMap<String, usize> = JavaHashMap::new();
        let assemblies: Vec<AsmP>;
        let primary_assemblies: Vec<AsmP>;
        if let Some(phase_set) = aa.phase_set() {
            let mut all = phase_set.all_assemblies();
            primary_assemblies = aa.assemblies();
            // Java (A39): only the sets merged straight into this one are counted, not those merged into
            // them before. Corrected: every level (PhaseSetMerger moves all their reads into this frame).
            let merged_phase_sets = if compat::java("A39") { phase_set.merged_phase_sets().clone() } else { phase_set.merged_phase_sets_deep() };
            for merged_phase_set in merged_phase_sets {
                all.extend(merged_phase_set.all_assemblies());
            }
            assemblies = all;
        } else {
            assemblies = aa.assemblies();
            primary_assemblies = assemblies.clone();
        }
        for assembly in &assemblies {
            let in_primary_assembly = primary_assemblies.contains(assembly);
            for read in assembly.support().clone() {
                let index = match fragment_map.get(&read.id().to_string()) {
                    Some(i) => *i,
                    None => {
                        fragment_map.insert(read.id().to_string(), fragments.len());
                        fragments.push(FragmentReads::new(read, in_primary_assembly));
                        continue;
                    }
                };
                let fragment_reads = &mut fragments[index];
                fragment_reads.in_primary_assembly |= in_primary_assembly;
                if fragment_reads.processed {
                    fragment_reads.set_read_details(&read);
                    continue;
                }
                fragment_reads.add_read(read);
                let first_read = fragment_reads.find_specific_read(true, false);
                let second_read = fragment_reads.find_specific_read(false, false);
                if first_read.is_none() || second_read.is_none() {
                    continue;
                }
                fragment_reads.processed = true;
                self.process_support_reads(first_read.as_ref(), second_read.as_ref(), in_primary_assembly);
                fragment_reads.set_fragment_length();
            }
        }
        let order: Vec<usize> = fragment_map.values().copied().collect();
        for index in order {
            let fragment_reads = &mut fragments[index];
            self.update_unique_fragment_positions(fragment_reads);
            if fragment_reads.processed {
                continue;
            }
            let first_read = fragment_reads.find_specific_read(true, true);
            let second_read = fragment_reads.find_specific_read(false, true);
            self.process_support_reads(first_read.as_ref(), second_read.as_ref(), fragment_reads.in_primary_assembly);
            fragment_reads.set_fragment_length();
        }
    }

    fn process_support_reads(&self, first_read: Option<&SupP>, second_read: Option<&SupP>, in_primary_assembly: bool) {
        let mut first_breakend_matches = first_read.map(|r| self.find_read_breakend_match(r)).unwrap_or_default();
        let mut second_breakend_matches = second_read.map(|r| self.find_read_breakend_match(r)).unwrap_or_default();
        if let (Some(first), Some(second)) = (first_read, second_read) {
            if !first_breakend_matches.is_empty() && second_breakend_matches.is_empty() && can_use_other_read_breakend_matches(first, second) {
                second_breakend_matches = first_breakend_matches.clone();
            } else if first_breakend_matches.is_empty() && !second_breakend_matches.is_empty() && can_use_other_read_breakend_matches(second, first) {
                first_breakend_matches = second_breakend_matches.clone();
            }
        }
        if !first_breakend_matches.is_empty() && !second_breakend_matches.is_empty() {
            self.process_complete_fragment(first_read.unwrap(), second_read.unwrap(), in_primary_assembly, &first_breakend_matches, &second_breakend_matches);
        } else if !first_breakend_matches.is_empty() {
            self.process_solo_read(first_read.unwrap(), in_primary_assembly, &first_breakend_matches);
        } else if !second_breakend_matches.is_empty() {
            self.process_solo_read(second_read.unwrap(), in_primary_assembly, &second_breakend_matches);
        }
    }

    fn process_complete_fragment(&self, first_read: &SupP, second_read: &SupP, in_primary_assembly: bool, first_breakend_matches: &[ReadBreakendMatch], second_breakend_matches: &[ReadBreakendMatch]) {
        let lower_index = first_read.full_assembly_index_start().min(second_read.full_assembly_index_start());
        let upper_index = first_read.full_assembly_index_end().max(second_read.full_assembly_index_end());
        let fragment_length = upper_index - lower_index + 1;
        first_read.m().set_inferred_fragment_length(fragment_length);
        second_read.m().set_inferred_fragment_length(fragment_length);
        let mut breakends: IdHashSet<BreakendP> = IdHashSet::new();
        add_unique_breakends(&mut breakends, first_breakend_matches);
        add_unique_breakends(&mut breakends, second_breakend_matches);
        let is_short_indel = breakends.iter().all(|x| x.is_short_local_del_dup_ins());
        let is_split_support = first_breakend_matches.iter().any(|x| x.is_split) || second_breakend_matches.iter().any(|x| x.is_split);
        if !is_split_support && is_short_indel {
            return;
        }
        let support_type = if is_split_support { SupportType::Junction } else { SupportType::Discordant };
        first_read.m().set_breakend_support_type(Some(support_type));
        second_read.m().set_breakend_support_type(Some(support_type));
        let sample_index = first_read.sample_index();
        for breakend in breakends.to_vec() {
            let mut targets = vec![breakend];
            if let Some(other) = breakend.other_breakend() {
                targets.push(other);
            }
            for b in targets {
                b.m().update_breakend_support(sample_index, is_split_support);
                b.m().update_breakend_read_strand(sample_index, first_read);
                b.m().update_breakend_read_strand(sample_index, second_read);
                b.m().add_inferred_fragment_length(fragment_length, true);
                if !in_primary_assembly {
                    b.m().add_non_primary_assembly_fragment_count();
                }
            }
        }
    }

    fn update_unique_fragment_positions(&self, fragment_reads: &FragmentReads) {
        if !fragment_reads.in_primary_assembly {
            return;
        }
        let mut processed_breakends: Vec<BreakendP> = Vec::new();
        for read in &fragment_reads.primary_reads {
            if !read.type_().is_split_support() {
                continue;
            }
            let breakend_matches = self.find_read_breakend_match(read);
            let fragment_position = read.five_prime_fragment_position();
            for breakend_match in breakend_matches {
                let breakend = breakend_match.breakend;
                if processed_breakends.contains(&breakend) {
                    continue;
                }
                if breakend.chromosome != read.chromosome() {
                    continue;
                }
                if !position_within(breakend.position, read.unclipped_start(), read.unclipped_end()) {
                    continue;
                }
                breakend.m().add_fragment_position(fragment_position);
                breakend.m().add_orientation_positions(read.orientation(), read.untrimmed_start(), read.untrimmed_end());
                processed_breakends.push(breakend);
            }
        }
    }

    fn process_solo_read(&self, read: &SupP, in_primary_assembly: bool, read_breakend_matches: &[ReadBreakendMatch]) {
        let aa = self.assembly_alignment;
        let mut breakends: IdHashSet<BreakendP> = IdHashSet::new();
        add_unique_breakends(&mut breakends, read_breakend_matches);
        let mut inferred_fragment_length = -1;
        let is_short_indel = breakends.iter().all(|x| x.is_short_local_del_dup_ins());
        let mut sv_type: Option<StructuralVariantType> = None;
        let mut indel_length = 0;
        if self.is_local_indel() {
            let link = aa.phase_set().unwrap().assembly_links()[0];
            indel_length = link.length();
            sv_type = Some(link.sv_type());
        } else if !read.is_discordant() && breakends.len() <= 2 {
            let first = *breakends.iter().next().expect("breakend");
            indel_length = first.sv_length();
            sv_type = Some(first.sv_type());
        }
        if sv_type.is_some_and(is_indel) && indel_length != 0 && read.insert_size() > 0 {
            inferred_fragment_length = self.calc_indel_solo_read_fragment_length(read, sv_type.unwrap(), indel_length);
        } else if !read.is_paired_read() {
            inferred_fragment_length = read.insert_size().abs();
        }
        if aa.breakends().len() == 2 && aa.assemblies().len() == 1 {
            inferred_fragment_length = self.calc_discordant_solo_read_fragment_length(read, aa.breakends());
        }
        let set_valid_fragment_length = inferred_fragment_length <= max_observed_concordant_frag_length();
        if set_valid_fragment_length {
            read.m().set_inferred_fragment_length(inferred_fragment_length);
        }
        let is_split_support = read_breakend_matches.iter().any(|x| x.is_split);
        if !is_split_support && is_short_indel {
            return;
        }
        read.m().set_breakend_support_type(Some(if is_split_support { SupportType::Junction } else { SupportType::Discordant }));
        let sample_index = read.sample_index();
        for breakend in breakends.to_vec() {
            let mut targets = vec![breakend];
            if let Some(other) = breakend.other_breakend() {
                targets.push(other);
            }
            for b in targets {
                b.m().update_breakend_support(sample_index, is_split_support);
                b.m().update_breakend_read_strand(sample_index, read);
                b.m().add_inferred_fragment_length(inferred_fragment_length, set_valid_fragment_length);
                if !in_primary_assembly {
                    b.m().add_non_primary_assembly_fragment_count();
                }
            }
        }
    }

    fn calc_indel_solo_read_fragment_length(&self, read: &SupP, sv_type: StructuralVariantType, indel_length: i32) -> i32 {
        let aa = self.assembly_alignment;
        let mut fragment_length = read.insert_size().abs();
        if aa.is_local_indel_assembly() {
            let assemblies = aa.assemblies();
            if read.is_paired_read() && sv_type == StructuralVariantType::Del {
                let Some(ref_assembly) = assemblies.iter().find(|x| x.support_count() == 0) else { return 0 };
                if (ref_assembly.is_forward_junction() && read.orientation().is_reverse()) || (ref_assembly.is_reverse_junction() && read.orientation().is_forward()) {
                    fragment_length -= indel_length;
                    return fragment_length;
                }
            }
            // a null pointer exception in Java if no assembly has support
            let junc_assembly = assemblies.iter().find(|x| x.support_count() > 0).expect("assembly with support");
            if junc_assembly.is_forward_junction() {
                fragment_length += read.right_clip_length();
            } else {
                fragment_length += read.left_clip_length();
            }
        } else if read.is_paired_read() {
            if sv_type == StructuralVariantType::Del {
                fragment_length -= indel_length;
            } else {
                fragment_length += indel_length;
            }
        }
        fragment_length
    }

    fn is_local_indel(&self) -> bool {
        let aa = self.assembly_alignment;
        let assemblies = aa.assemblies();
        if assemblies.len() != 2 {
            return false;
        }
        let Some(phase_set) = aa.phase_set() else { return false };
        if phase_set.assembly_links().len() != 1 {
            return false;
        }
        if assemblies.iter().all(|x| x.indel()) {
            return true;
        }
        if assemblies.iter().all(|x| x.outcome() == AssemblyOutcome::LocalIndel) {
            return true;
        }
        let link = phase_set.assembly_links()[0];
        is_indel(link.sv_type()) && link.length() <= 1000
    }

    fn calc_discordant_solo_read_fragment_length(&self, read: &SupP, breakends: &[BreakendP]) -> i32 {
        let original_assembly = self.assembly_alignment.assemblies()[0];
        let mut local_breakend: Option<BreakendP> = None;
        let mut remote_breakend: Option<BreakendP> = None;
        let mut remote_fragment_length_part = -1;
        for breakend in breakends {
            if breakend.orient == original_assembly.junction().orient && breakend.chromosome == read.chromosome() && position_within(breakend.position, read.unclipped_start(), read.unclipped_end()) {
                local_breakend = Some(*breakend);
                continue;
            }
            if Some(breakend.chromosome.as_str()) == read.mate_chromosome() && breakend.orient == read.mate_orientation() {
                if breakend.orient.is_forward() && read.mate_alignment_end() <= breakend.position {
                    remote_breakend = Some(*breakend);
                    remote_fragment_length_part = breakend.position - read.mate_alignment_start();
                } else if breakend.orient.is_reverse() && read.mate_alignment_start() >= breakend.position {
                    remote_breakend = Some(*breakend);
                    remote_fragment_length_part = read.mate_alignment_end() - breakend.position;
                }
            }
        }
        if let (Some(_), Some(local)) = (remote_breakend, local_breakend) {
            let adjustment = local.inserted_bases.len() as i32 - local.homology.length();
            return read.junction_read_start_distance() + remote_fragment_length_part + adjustment;
        }
        if !read.is_discordant() && read.type_() == SupportType::Junction {
            return read.insert_size() + read.left_clip_length().max(read.right_clip_length());
        }
        -1
    }

    fn find_read_breakend_match(&self, read: &SupP) -> Vec<ReadBreakendMatch> {
        let mut closest_discordant_breakend: Option<BreakendP> = None;
        let mut closest_disc_junction_distance = INVALID_DISCORANT_DISTANCE;
        let mut breakend_matches: Vec<ReadBreakendMatch> = Vec::with_capacity(2);
        for breakend in self.assembly_alignment.breakends() {
            if self.read_spans_junction(breakend, read) {
                breakend_matches.push(ReadBreakendMatch { breakend: *breakend, is_split: true });
                continue;
            }
            let discordant_distance = read_discordant_breakend_distance(breakend, read);
            if discordant_distance != INVALID_DISCORANT_DISTANCE && (closest_disc_junction_distance == INVALID_DISCORANT_DISTANCE || discordant_distance < closest_disc_junction_distance) {
                closest_disc_junction_distance = discordant_distance;
                closest_discordant_breakend = Some(*breakend);
            }
        }
        if let Some(b) = closest_discordant_breakend {
            breakend_matches.push(ReadBreakendMatch { breakend: b, is_split: false });
        }
        breakend_matches
    }

    fn read_spans_junction(&self, breakend: &Breakend, read: &SupportRead) -> bool {
        if read_alignment_spans_junction(breakend, read) {
            return true;
        }
        let read_seq_index_start = read.full_assembly_index_start();
        let read_seq_index_end = read.full_assembly_index_end();
        let full_sequence_end_index = self.assembly_alignment.full_sequence_length() - 1;
        for segment in breakend.segments() {
            if let Some(indel_sequence_indices) = segment.indel_sequence_indices() {
                if read_seq_index_start < indel_sequence_indices[0] && read_seq_index_end >= indel_sequence_indices[0] {
                    return true;
                }
                if read_seq_index_start < indel_sequence_indices[1] && read_seq_index_end >= indel_sequence_indices[1] {
                    return true;
                }
            }
            let segment_junction_index = if breakend.orient.is_forward() {
                if segment.alignment.is_forward() { segment.alignment.sequence_end() } else { segment.alignment.sequence_start() }
            } else if segment.alignment.is_forward() {
                segment.alignment.sequence_start()
            } else {
                segment.alignment.sequence_end()
            };
            if segment_junction_index == 0 || segment_junction_index == full_sequence_end_index {
                continue;
            }
            if read_seq_index_start < segment_junction_index && read_seq_index_end > segment_junction_index {
                return true;
            }
        }
        false
    }
}
