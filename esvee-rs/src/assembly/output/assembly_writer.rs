//! AssemblyWriter and AssemblyWriterUtils: the assembly TSV.

use std::collections::HashSet;
use std::io::Write;

use crate::assembly::assembly_utils::has_unset_bases;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::globals::READ_ID_TRIMMER;
use crate::assembly::types::assembly_link::LinkP;
use crate::assembly::types::assembly_stats::AssemblyStats;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::ref_side_soft_clip::RefSideSoftClipP;
use crate::assembly::types::remote_region::RemoteRegionP;
use crate::assembly::types::repeat_info::RepeatInfo;
use crate::assembly::types::{LinkType, RemoteReadType, SupportType};
use crate::common::java_sort::java_sort;

pub struct AssemblyWriter {
    writer: Option<std::io::BufWriter<std::fs::File>>,
    assembly_detailed_tsv: bool,
    assembly_ref_base_write_max: i32,
    saga: bool,
}

impl AssemblyWriter {
    pub fn new(config: &AssemblyConfig) -> AssemblyWriter {
        let mut w = AssemblyWriter { writer: None, assembly_detailed_tsv: config.assembly_detailed_tsv, assembly_ref_base_write_max: config.assembly_ref_base_write_max, saga: config.saga_fasta_file.is_some() };
        if !config.write_types.contains(&WriteType::JuncAssembly) {
            return w;
        }
        match std::fs::File::create(config.output_filename(WriteType::JuncAssembly)) {
            Ok(f) => {
                let mut out = std::io::BufWriter::new(f);
                let mut sj: Vec<String> = Vec::new();
                for h in ["Id", "Chromosome", "JuncPosition", "JuncOrientation", "JuncType", "ExtBaseLength", "RefBasePosition", "RefBaseLength", "RefBaseCigar"] {
                    sj.push(h.to_string());
                }
                add_support_header(&mut sj);
                sj.push("Outcome".to_string());
                add_phasing_header(&mut sj);
                for h in ["JunctionSequence", "RefBaseSequence", "InsertType", "RefBaseCandidates", "UnmappedCandidates", "AssemblyInfo", "RawDiscPosition"] {
                    sj.push(h.to_string());
                }
                if w.assembly_detailed_tsv {
                    for h in ["InitialReadId", "ExtBaseBuildInfo", "MismatchReads", "RefSideSoftClips", "RefBaseTrimmed", "RefBaseTrimLength", "RepeatInfo"] {
                        sj.push(h.to_string());
                    }
                    AssemblyStats::add_read_stats_header(&mut sj);
                    AssemblyStats::add_read_type_header(&mut sj);
                    add_remote_region_header(&mut sj);
                    sj.push("MergedAssemblies".to_string());
                }
                if w.saga {
                    sj.push("SagaMatch".to_string());
                    sj.push("SagaMatchCigar".to_string());
                }
                let _ = writeln!(out, "{}", sj.join("\t"));
                w.writer = Some(out);
            }
            Err(e) => crate::error!("failed to initialise assembly writer: {}", e),
        }
        w
    }

    pub fn enabled(&self) -> bool {
        self.writer.is_some()
    }

    pub fn close(&mut self) {
        if let Some(w) = self.writer.as_mut() {
            let _ = w.flush();
        }
    }

    pub fn write_assembly(&mut self, assembly: AsmP) {
        if self.writer.is_none() {
            return;
        }
        let junction = assembly.junction();
        let mut sj: Vec<String> = Vec::new();
        sj.push(assembly.id().to_string());
        sj.push(junction.chromosome.to_string());
        sj.push(junction.position.to_string());
        sj.push(junction.orient.to_string());
        sj.push((if junction.discordant_only { "DISC" } else if junction.indel_based() { "INDEL" } else { "SPLIT" }).to_string());
        sj.push(assembly.extension_length().to_string());
        sj.push(assembly.ref_base_position().to_string());
        sj.push(assembly.ref_base_length().to_string());
        sj.push(assembly.ref_base_cigar());
        add_support_counts(&assembly, &mut sj);
        sj.push(assembly.outcome().name().to_string());
        add_phasing_info(assembly, &mut sj);
        if has_unset_bases(&assembly) {
            sj.push("UNSET_BASES".to_string());
            sj.push("UNSET_BASES".to_string());
        } else {
            sj.push(assembly.form_junction_sequence());
            let ref_base_length = if self.assembly_ref_base_write_max == 0 { assembly.ref_base_length() } else { self.assembly_ref_base_write_max };
            sj.push(assembly.form_ref_base_sequence_max(ref_base_length));
        }
        sj.push((if assembly.has_line_sequence() { "LINE" } else { "NONE" }).to_string());
        sj.push(assembly.stats().candidate_support_count.to_string());
        sj.push(assembly.stats().unmapped_read_count.to_string());
        sj.push(assembly.assembly_alignment_info());
        let raw_discordant_position = junction.raw_discordant_position();
        sj.push(if raw_discordant_position > 0 { raw_discordant_position.to_string() } else { String::new() });
        if self.assembly_detailed_tsv {
            sj.push(match assembly.initial_read_id() {
                Some(id) => READ_ID_TRIMMER.restore(id),
                None => {
                    if READ_ID_TRIMMER.enabled() {
                        READ_ID_TRIMMER.restore("null")
                    } else {
                        "null".to_string()
                    }
                }
            });
            sj.push(assembly.ext_base_build_info().to_string());
            sj.push(assembly.mismatch_read_count().to_string());
            sj.push(ref_side_soft_clips_str(assembly.ref_side_soft_clips()));
            sj.push(assembly.ref_bases_repeated_trimmed().unwrap_or("null").to_string());
            sj.push(assembly.ref_base_trim_length().to_string());
            sj.push(repeats_info_str(assembly.repeat_info()));
            assembly.stats().add_read_stats(&mut sj);
            assembly.stats().add_read_type_counts(&mut sj);
            add_remote_region_info(assembly, &mut sj);
            sj.push(assembly.merged_assembly_count().to_string());
        }
        if self.saga {
            match assembly.saga_match() {
                Some(m) => {
                    sj.push(m.variant_string());
                    sj.push(m.cigar_string());
                }
                None => {
                    sj.push(String::new());
                    sj.push(String::new());
                }
            }
        }
        let _ = writeln!(self.writer.as_mut().unwrap(), "{}", sj.join("\t"));
    }
}

pub fn add_support_header(sj: &mut Vec<String>) {
    for h in ["SplitFrags", "RefSplitFrags", "DiscFrags", "RefDiscFrags"] {
        sj.push(h.to_string());
    }
}

pub fn add_support_counts(assembly: &JunctionAssembly, sj: &mut Vec<String>) {
    let mut junc_frags = 0;
    let mut ref_sample_junc_frags = 0;
    let mut disc_frags = 0;
    let mut ref_sample_disc_frags = 0;
    let mut unique_read_ids: HashSet<&str> = HashSet::new();
    for read in assembly.support() {
        let is_reference = read.is_reference();
        if unique_read_ids.contains(read.id()) {
            continue;
        }
        if read.type_() != SupportType::Extension {
            unique_read_ids.insert(read.id());
        }
        if read.type_().is_split_support() {
            junc_frags += 1;
            if is_reference {
                ref_sample_junc_frags += 1;
            }
        } else if read.type_() == SupportType::Discordant {
            disc_frags += 1;
            if is_reference {
                ref_sample_disc_frags += 1;
            }
        }
    }
    for v in [junc_frags, ref_sample_junc_frags, disc_frags, ref_sample_disc_frags] {
        sj.push(v.to_string());
    }
}

pub fn add_phasing_header(sj: &mut Vec<String>) {
    for h in ["PhaseGroupId", "PhaseGroupCount", "PhaseSetId", "PhaseSetCount", "SplitLinks", "FacingLinks", "SvType", "SvLength", "InsertedBases", "OverlapBases", "SecondaryLinks"] {
        sj.push(h.to_string());
    }
}

pub fn add_phasing_info(assembly: AsmP, sj: &mut Vec<String>) {
    let phase_group = assembly.phase_group();
    let phase_group_id = phase_group.map(|pg| pg.id()).unwrap_or(-1);
    let phase_group_count = phase_group.map(|pg| pg.assembly_count()).unwrap_or(1);
    let phase_set_id;
    let phase_set_count;
    let mut split_link_info = String::new();
    let mut facing_link_info = String::new();
    let mut sv_type = String::new();
    let mut sv_link_length = 0;
    let mut inserted_bases = String::new();
    let mut overlap_bases = String::new();
    let mut secondary_link_info = String::new();
    match JunctionAssembly::phase_set(assembly) {
        Some(phase_set) => {
            phase_set_id = if phase_set.merged() { format!("{}_{}", phase_set.id(), phase_set.merged_phase_set_id()) } else { phase_set.id().to_string() };
            phase_set_count = phase_set.assemblies().len();
            let assembly_links = phase_set.find_assembly_links(assembly);
            let split_links: Vec<LinkP> = assembly_links.iter().filter(|x| x.type_() == LinkType::Split || x.type_() == LinkType::Indel).copied().collect();
            let facing_links: Vec<LinkP> = assembly_links.iter().filter(|x| x.type_() == LinkType::Facing).copied().collect();
            split_link_info = assembly_links_str(assembly, &split_links);
            facing_link_info = assembly_links_str(assembly, &facing_links);
            if let Some(sv_link) = split_links.first() {
                sv_type = sv_link.sv_type().name().to_string();
                sv_link_length = sv_link.length();
                inserted_bases = sv_link.inserted_bases().to_string();
                overlap_bases = sv_link.overlap_bases().to_string();
            }
            let secondary_split_links = phase_group.expect("phase group").find_secondary_split_links(assembly);
            secondary_link_info = assembly_links_str(assembly, &secondary_split_links);
        }
        None => {
            phase_set_id = assembly.non_phase_set_id().to_string();
            phase_set_count = 1;
        }
    }
    sj.push(phase_group_id.to_string());
    sj.push(phase_group_count.to_string());
    sj.push(phase_set_id);
    sj.push(phase_set_count.to_string());
    sj.push(split_link_info);
    sj.push(facing_link_info);
    sj.push(sv_type);
    sj.push(sv_link_length.to_string());
    sj.push(inserted_bases);
    sj.push(overlap_bases);
    sj.push(secondary_link_info);
}

pub fn assembly_links_str(assembly: AsmP, assembly_links: &[LinkP]) -> String {
    assembly_links.iter().map(|link| link.other_assembly(assembly).junction().coords()).collect::<Vec<_>>().join(";")
}

pub fn repeats_info_str(repeats: &[RepeatInfo]) -> String {
    if repeats.is_empty() {
        return String::new();
    }
    let mut longest: Option<&RepeatInfo> = None;
    let mut longest_sequence: Option<&RepeatInfo> = None;
    for repeat in repeats {
        if longest.is_none() || repeat.count > longest.unwrap().count {
            longest = Some(repeat);
        }
        // as in Java, compared with the repeat of the highest count
        if longest_sequence.is_none() || repeat.total_length() > longest.unwrap().total_length() {
            longest_sequence = Some(repeat);
        }
    }
    let (l, ls) = (longest.unwrap(), longest_sequence.unwrap());
    format!("{} max({}={}) long({}={})", repeats.len(), l.bases, l.count, ls.bases, ls.count)
}

pub fn add_remote_region_header(sj: &mut Vec<String>) {
    for h in ["RemoteRegionCount", "RemoteRegionInfo", "RemoteRegionJuncMate", "RemoteRegionJuncSupps", "RemoteRegionDiscordant"] {
        sj.push(h.to_string());
    }
}

pub fn add_remote_region_info(assembly: AsmP, sj: &mut Vec<String>) {
    sj.push(assembly.remote_regions().len().to_string());
    sj.push(remote_region_info_str(assembly.m().remote_regions_mut()));
    let mut remote_junct_supp = 0;
    let mut remote_junct_mate = 0;
    let mut remote_discordant = 0;
    for region in assembly.remote_regions() {
        remote_junct_mate += region.read_type_counts()[RemoteReadType::JunctionMate.ordinal()];
        remote_junct_supp += region.read_type_counts()[RemoteReadType::Supplementary.ordinal()];
        remote_discordant += region.read_type_counts()[RemoteReadType::Discordant.ordinal()];
    }
    sj.push(remote_junct_mate.to_string());
    sj.push(remote_junct_supp.to_string());
    sj.push(remote_discordant.to_string());
}

/// sorts the assembly's list of regions by read count, as in Java
fn remote_region_info_str(regions: &mut Vec<RemoteRegionP>) -> String {
    if regions.is_empty() {
        return String::new();
    }
    java_sort(regions, |a, b| (-a.read_count()).cmp(&(-b.read_count())));
    regions.iter().take(3).map(|region| format!("{}:{}-{}={}", region.chromosome, region.start(), region.end(), region.read_count())).collect::<Vec<_>>().join(" ")
}

pub fn ref_side_soft_clips_str(ref_side_soft_clips: &[RefSideSoftClipP]) -> String {
    ref_side_soft_clips.iter().map(|x| format!("{}:{}={}", x.position, x.max_length(), x.read_count())).collect::<Vec<_>>().join(";")
}
