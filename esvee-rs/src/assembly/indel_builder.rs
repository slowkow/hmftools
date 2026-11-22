//! IndelBuilder: indel-based junctions and reads whose indels are treated as soft clips.

use std::collections::HashMap;

use crate::assembly::common_utils::IndelCoords;
use crate::assembly::constants::*;
use crate::assembly::ptr::A;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::{get_read_index_at_reference_position, record_soft_clips_at_junction};
use crate::assembly::sam::CigarOp;
use crate::assembly::types::junction::{JuncP, Junction};
use crate::assembly::types::junction_assembly::JunctionAssembly;
use crate::assembly::types::{AssemblyOutcome, SupportType};
use crate::common::java_hash::JavaHashMap;
use crate::common::java_sort::java_sort;

pub fn set_indel_status(junction: JuncP, reads: &[ReadP]) {
    let mut soft_clipped_count = 0;
    let mut matched_indel_count = 0;
    for read in reads {
        let ic = read.indel_coords();
        if ic.as_ref().is_some_and(|c| c.length >= MIN_INDEL_LENGTH && c.matches_junction(junction.position, junction.orient)) {
            matched_indel_count += 1;
        } else if record_soft_clips_at_junction(read, &junction) {
            // as in Java, both tests are on a forward junction
            if junction.is_forward() && read.right_clip_length() >= MIN_VARIANT_LENGTH {
                soft_clipped_count += 1;
            } else if junction.is_forward() && read.left_clip_length() >= MIN_VARIANT_LENGTH {
                soft_clipped_count += 1;
            }
        }
    }
    if junction.indel_based() {
        if matched_indel_count < ASSEMBLY_MIN_READ_SUPPORT && soft_clipped_count > ASSEMBLY_MIN_READ_SUPPORT {
            junction.m().set_indel_status(false);
        } else if soft_clipped_count > matched_indel_count * 2 {
            junction.m().set_indel_status(false);
        }
    } else if soft_clipped_count < ASSEMBLY_MIN_READ_SUPPORT && matched_indel_count > ASSEMBLY_MIN_READ_SUPPORT {
        junction.m().set_indel_status(true);
    } else if matched_indel_count > soft_clipped_count * 2 {
        junction.m().set_indel_status(true);
    }
}

pub fn find_indel_extension_reads(junction: JuncP, raw_reads: &[ReadP], extension_reads: &mut Vec<ReadP>, junction_reads: &mut Vec<ReadP>, non_junction_reads: &mut Vec<ReadP>) {
    let mut indel_coords_count: HashMap<i32, i32> = HashMap::new();
    let mut dominant_indel: Option<A<IndelCoords>> = None;
    let mut dominant_indel_freq = 0;
    for read in raw_reads {
        if read.invalid_indel() {
            continue;
        }
        if let Some(indel_coords) = read.indel_coords() {
            if !indel_coords.matches_junction(junction.position, junction.orient) {
                continue;
            }
            if indel_coords.length >= MIN_INDEL_LENGTH {
                extension_reads.push(read.clone());
                let freq = indel_coords_count.entry(indel_coords.length).or_insert(0);
                *freq += 1;
                if *freq > dominant_indel_freq {
                    dominant_indel_freq = *freq;
                    dominant_indel = Some(indel_coords);
                }
            } else {
                junction_reads.push(read.clone());
            }
        } else {
            let mut crosses_junction = false;
            if junction.is_forward() {
                if read.has_indel_implied_unclipped_end() && read.max_unclipped_end() > junction.position {
                    crosses_junction = true;
                } else if read.is_right_clipped() && read.unclipped_end() > junction.position {
                    crosses_junction = true;
                }
            } else if read.has_indel_implied_unclipped_start() && read.min_unclipped_start() < junction.position {
                crosses_junction = true;
            } else if read.is_left_clipped() && read.unclipped_start() < junction.position {
                crosses_junction = true;
            }
            if crosses_junction {
                junction_reads.push(read.clone());
            } else {
                non_junction_reads.push(read.clone());
            }
        }
    }
    if dominant_indel.is_some() {
        junction.m().set_indel_coords(dominant_indel);
    }
}

pub fn calc_indel_inferred_unclipped_positions(read: &mut Read) -> bool {
    let n = read.cigar_elements().len();
    if n < 3 {
        return false;
    }
    let mut insert_total = 0;
    let mut aligned_total = 0;
    let mut has_min_length_indel = false;
    let last_index = n - 1;
    for i in 0..n {
        let element = read.cigar_elements()[i];
        if element.op.is_indel() {
            has_min_length_indel |= element.len >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
            if i == 0 || (i == 1 && read.cigar_elements()[0].op != CigarOp::M) {
                read.mark_invalid_indel();
                return false;
            }
            if i == last_index || (i == last_index - 1 && read.cigar_elements()[last_index].op != CigarOp::M) {
                read.mark_invalid_indel();
                return false;
            }
            if element.op == CigarOp::I {
                insert_total += element.len;
            }
        } else if element.op == CigarOp::M {
            aligned_total += element.len;
        }
    }
    if !has_min_length_indel {
        return false;
    }
    let mut unclipped_start = None;
    let mut unclipped_end = None;
    if !read.is_left_clipped() {
        unclipped_start = Some(read.alignment_end() - aligned_total + 1 - insert_total);
    }
    if !read.is_right_clipped() {
        unclipped_end = Some(read.alignment_start() + aligned_total - 1 + insert_total);
    }
    read.set_indel_inferred_unclipped_positions(unclipped_start, unclipped_end);
    true
}

pub fn converted_indel_crosses_junction(junction: &Junction, read: &Read) -> bool {
    if junction.is_forward() {
        if read.unclipped_end() > junction.position {
            return read.has_indel_implied_unclipped_end();
        }
    } else if read.unclipped_start() < junction.position {
        return read.has_indel_implied_unclipped_start();
    }
    false
}

pub fn indel_assemblies_match(first: &JunctionAssembly, second: &JunctionAssembly) -> bool {
    match (first.indel_coords(), second.indel_coords()) {
        (Some(a), Some(b)) if first.is_forward_junction() != second.is_forward_junction() => a.matches(b),
        _ => false,
    }
}

pub fn find_inserted_bases(read: &Read) -> String {
    let indel_coords = read.indel_coords().expect("read without indel coords");
    let mut read_index = 0usize;
    let mut ref_position = read.alignment_start();
    for element in read.cigar_elements() {
        if (element.op == CigarOp::I && element.len == indel_coords.length) || ref_position == indel_coords.pos_end {
            let mut inserted_bases = String::new();
            for i in 0..element.len as usize {
                inserted_bases.push(read.get_bases()[read_index + i] as char);
            }
            return inserted_bases;
        }
        if element.op.consumes_read_bases() {
            read_index += element.len as usize;
        }
        if element.op.consumes_reference_bases() {
            ref_position += element.len;
        }
    }
    String::new()
}

pub fn check_matches_local_indel_reads(reads: &[ReadP], junction: &Junction, extension_bases: &[u8]) -> i32 {
    let mut indel_reads_map: JavaHashMap<String, Vec<ReadP>> = JavaHashMap::new();
    for read in reads {
        let Some(indel_coords) = read.indel_coords() else { continue };
        if read.left_clip_length() > 0 || read.right_clip_length() > 0 {
            continue;
        }
        if (junction.is_forward() && indel_coords.pos_end > junction.position) || (junction.is_reverse() && indel_coords.pos_start < junction.position) {
            continue;
        }
        let key = indel_coords.to_string();
        if !indel_reads_map.contains_key(&key) {
            indel_reads_map.insert(key.clone(), Vec::new());
        }
        indel_reads_map.get_mut(&key).unwrap().push(read.clone());
    }
    if indel_reads_map.is_empty() {
        return 0;
    }
    let mut top_indel_coords_str: Option<String> = None;
    let mut top_frequency = 0;
    for (key, value) in indel_reads_map.iter() {
        if value.len() > top_frequency {
            top_frequency = value.len();
            top_indel_coords_str = Some(key.clone());
        }
    }
    let total_indel_reads: usize = indel_reads_map.iter_unordered().map(|(_, v)| v.len()).sum();
    if (top_frequency as f64) < total_indel_reads as f64 * 0.5 {
        return 0;
    }
    let indel_reads = indel_reads_map.get(&top_indel_coords_str.unwrap()).unwrap();
    let mut matching_read_count = 0;
    for read in indel_reads {
        if read_has_extension_sequence(read, junction, extension_bases) {
            matching_read_count += 1;
        }
    }
    matching_read_count
}

fn read_has_extension_sequence(read: &Read, junction: &Junction, extension_bases: &[u8]) -> bool {
    let indel_coords = read.indel_coords().expect("read without indel coords");
    let is_forward = junction.is_forward();
    let indel_position = if is_forward { indel_coords.pos_start } else { indel_coords.pos_end };
    let indel_read_index = get_read_index_at_reference_position(read, indel_position, false);
    let indel_junction_distance = (indel_position - junction.position).abs();
    let mut implied_junction_read_index;
    if is_forward {
        implied_junction_read_index = indel_read_index + indel_junction_distance;
        if indel_coords.is_insert() {
            implied_junction_read_index -= indel_coords.length;
        } else {
            implied_junction_read_index += indel_coords.length;
        }
    } else {
        implied_junction_read_index = indel_read_index - indel_junction_distance;
        if indel_coords.is_insert() {
            implied_junction_read_index += indel_coords.length;
        } else {
            implied_junction_read_index -= indel_coords.length;
        }
    }
    let mut overlap_bases = 0;
    let mut ext_seq_index: i32 = if is_forward { 0 } else { extension_bases.len() as i32 - 1 };
    let mut read_index = implied_junction_read_index;
    while read_index >= 0 && read_index < read.bases_length() {
        // an index outside the extension bases is an ArrayIndexOutOfBoundsException in Java
        if read.get_bases()[read_index as usize] != extension_bases[ext_seq_index as usize] {
            break;
        }
        overlap_bases += 1;
        if overlap_bases >= ASSEMBLY_READ_OVERLAP_BASES {
            return true;
        }
        if is_forward {
            ext_seq_index += 1;
            read_index += 1;
        } else {
            ext_seq_index -= 1;
            read_index -= 1;
        }
    }
    false
}

pub fn is_weak_indel_based_unlinked_assembly(assembly: &JunctionAssembly) -> bool {
    if assembly.junction().indel_based() {
        return false;
    }
    if assembly.outcome() != AssemblyOutcome::Unset && assembly.outcome() != AssemblyOutcome::NoLink {
        return false;
    }
    if assembly.extension_length() >= ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_MIN_LENGTH {
        return false;
    }
    let mut total_junc_reads = 0;
    let mut max_non_indel_length = 0;
    let mut indel_read_lengths: Vec<i32> = Vec::new();
    for read in assembly.support() {
        if read.type_().is_split_support() {
            let extension_length = read.extension_base_matches();
            total_junc_reads += 1;
            if read.type_() == SupportType::Indel {
                indel_read_lengths.push(extension_length);
            } else {
                max_non_indel_length = max_non_indel_length.max(extension_length);
            }
        }
    }
    if indel_read_lengths.len() as f64 >= ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_INDEL_PERC * total_junc_reads as f64 {
        return true;
    }
    if indel_read_lengths.len() < 2 {
        return false;
    }
    java_sort(&mut indel_read_lengths, |a, b| a.cmp(b));
    let second_longest_indel = indel_read_lengths[indel_read_lengths.len() - 2];
    second_longest_indel > max_non_indel_length
}
