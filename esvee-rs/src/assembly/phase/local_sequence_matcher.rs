//! LocalSequenceMatcher: looks for an assembly's extension bases in the reference close to its junction,
//! which makes a short local DEL or DUP.

use crate::assembly::assembly_utils::create_base_quals_above_min_threshold;
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::phase::assembly_linker::find_best_sequence_match;
use crate::assembly::ptr::P;
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::sequence_compare::index_of;
use crate::assembly::types::assembly_link::{AssemblyLink, LinkP};
use crate::assembly::types::junction::Junction;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::junction_sequence::JunctionSequence;
use crate::assembly::types::LinkType;

pub struct LocalSequenceMatcher<'a> {
    ref_genome: Option<&'a dyn RefGenomeInterface>,
    local_match_distance: i32,
}

fn index_of_from(haystack: &[u8], needle: &[u8], from: i32) -> i32 {
    let from = from.max(0) as usize;
    if from > haystack.len() {
        return -1;
    }
    index_of(&haystack[from..], needle).map(|i| i + from as i32).unwrap_or(-1)
}

impl<'a> LocalSequenceMatcher<'a> {
    pub fn new(ref_genome: Option<&'a dyn RefGenomeInterface>, local_match_distance: i32) -> Self {
        LocalSequenceMatcher { ref_genome, local_match_distance }
    }

    pub fn try_local_assembly_link(&self, assembly: AsmP) -> Option<LinkP> {
        let ref_genome = self.ref_genome.expect("ref genome");
        let junction = assembly.junction();
        let local_region_start = junction.position - self.local_match_distance;
        let local_region_end = junction.position + self.local_match_distance;
        let chromosome_end = ref_genome.get_chromosome_length(junction.chromosome);
        if local_region_start < 1 || local_region_end >= chromosome_end {
            return None;
        }
        let ref_genome_bases = ref_genome.get_bases(junction.chromosome, local_region_start, local_region_end).expect("ref genome bases");
        let ref_base_quals = create_base_quals_above_min_threshold(ref_genome_bases.len());
        // Java (A40): the sequence compared with the local reference is the straddling one, up to 50 ref
        // bases and then the whole extension, while the seed positions found below are offsets in the
        // extension alone; find_best_sequence_match takes them as offsets in that longer sequence. At a
        // forward junction the comparison is then shifted by the ref bases taken (up to 50), so with any
        // mismatch in the extension (no exact match below) a local link is only found where those ref
        // bases happen to repeat the local reference. At a reverse junction the frames agree, but the 100
        // bases compared take ref bases after an extension shorter than 100, which by construction do not
        // follow the extension's local match, so the link is found only for extensions of about 97 bases
        // or more. Corrected: the extension alone is compared (Java's commented-out alternative), and the
        // link is formed at the local position matched to the first extension base, as for an exact match.
        let java_a40 = crate::common::compat::java("A40");
        let assembly_seq = if java_a40 {
            JunctionSequence::form_straddling_match_sequence(&assembly)
        } else {
            JunctionSequence::form_full_extension_match_sequence(&assembly, false)
        };
        let local_ref_orientation = junction.orient.opposite();
        let local_ref_seq = JunctionSequence::from_bases(&ref_genome_bases, &ref_base_quals, local_ref_orientation, false);
        let ref_full = local_ref_seq.full_sequence.as_bytes();

        let assembly_ext_bases_str = assembly.form_junction_sequence();
        let assembly_ext_bases = assembly_ext_bases_str.as_bytes();
        let assembly_junction_seq_length = assembly_ext_bases.len() as i32;
        let assembly_seq_index_in_ref = index_of_from(ref_full, assembly_ext_bases, 0);
        if assembly_seq_index_in_ref >= 0 {
            let next_start_index = assembly_seq_index_in_ref + assembly_junction_seq_length;
            if index_of_from(ref_full, assembly_ext_bases, next_start_index) >= 0 {
                return None;
            }
            return Some(self.form_local_link(assembly, local_region_start, assembly_seq_index_in_ref));
        }
        let mut alternative_index_starts: Vec<[i32; 2]> = Vec::new();
        let subsequence_length = MATCH_SUBSEQUENCE_LENGTH;
        let sub_seq_iterations = assembly_junction_seq_length / subsequence_length;
        for i in 0..sub_seq_iterations {
            let match_seq_start_index = i * subsequence_length;
            let match_seq_end_index = match_seq_start_index + subsequence_length;
            if match_seq_end_index >= assembly_junction_seq_length {
                break;
            }
            let assembly_sub_sequence = &assembly_ext_bases[match_seq_start_index as usize..match_seq_end_index as usize];
            let mut ref_sub_seq_index = index_of_from(ref_full, assembly_sub_sequence, 0);
            if ref_sub_seq_index < 0 {
                continue;
            }
            alternative_index_starts.push([match_seq_start_index, ref_sub_seq_index]);
            ref_sub_seq_index = index_of_from(ref_full, assembly_sub_sequence, ref_sub_seq_index + subsequence_length);
            while ref_sub_seq_index >= 0 {
                alternative_index_starts.push([match_seq_start_index, ref_sub_seq_index]);
                ref_sub_seq_index = index_of_from(ref_full, assembly_sub_sequence, ref_sub_seq_index + subsequence_length);
            }
        }
        let min_overlap_length = assembly.extension_length().min(ASSEMBLY_LINK_OVERLAP_BASES);
        let top_match_indices = find_best_sequence_match(&assembly_seq, &local_ref_seq, min_overlap_length, &alternative_index_starts)?;
        let local_ref_index_start = if java_a40 {
            top_match_indices[1]
        } else {
            // the match pairs an assembly index with a local index; move to the first extension base
            let extension_start_index = if assembly.is_forward_junction() { assembly.junction_index() + 1 } else { 0 };
            top_match_indices[1] + extension_start_index - top_match_indices[0]
        };
        Some(self.form_local_link(assembly, local_region_start, local_ref_index_start))
    }

    fn form_local_link(&self, assembly: AsmP, local_region_start: i32, local_ref_index_start: i32) -> LinkP {
        let mut local_ref_junction_pos = local_region_start + local_ref_index_start;
        let local_ref_junction_index;
        let local_ref_orientation;
        let extension_base_length = assembly.extension_length();
        let fixed_assembly_ref_base_length = 50;
        let assembly_ref_bases = assembly.form_ref_base_sequence_max(fixed_assembly_ref_base_length);
        let assembly_ext_bases = assembly.form_junction_sequence();
        if assembly.is_forward_junction() {
            local_ref_orientation = Orientation::Reverse;
            local_ref_junction_index = assembly_ref_bases.len() as i32;
        } else {
            local_ref_orientation = Orientation::Forward;
            local_ref_junction_pos += extension_base_length - 1;
            local_ref_junction_index = assembly_ext_bases.len() as i32 - 1;
        }
        let local_ref_junction = P::new(Junction::new(assembly.junction().chromosome, local_ref_junction_pos, local_ref_orientation, false, false, false, 0, None));
        let local_assembly_sequence = if local_ref_orientation.is_forward() { format!("{assembly_ext_bases}{assembly_ref_bases}") } else { format!("{assembly_ref_bases}{assembly_ext_bases}") };
        let local_assembly_bases = local_assembly_sequence.into_bytes();
        let local_assembly_quals = create_base_quals_above_min_threshold(local_assembly_bases.len());
        let local_ref_assembly = JunctionAssembly::new(local_ref_junction, local_assembly_bases, local_assembly_quals, &[], Vec::new());
        local_ref_assembly.m().set_junction_index(local_ref_junction_index);
        if assembly.is_forward_junction() {
            AssemblyLink::new(assembly, local_ref_assembly, LinkType::Split, "", "")
        } else {
            AssemblyLink::new(local_ref_assembly, assembly, LinkType::Split, "", "")
        }
    }
}
