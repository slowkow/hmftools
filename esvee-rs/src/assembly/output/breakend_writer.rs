//! BreakendWriter: the breakend TSV.

use std::collections::HashMap;
use std::io::{BufWriter, Write};

use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend::BreakendP;
use crate::assembly::common_utils::within_line_proximity;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::line_utils::is_mobile_line_element;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::InsertionType;
use crate::common::fmt::java_fixed;

const OTHER_BREAKEND_COORDS: &str = "_OTHER";

pub struct BreakendWriter {
    writer: Option<BufWriter<std::fs::File>>,
    tumor_count: usize,
    saga: bool,
}

impl BreakendWriter {
    pub fn new(config: &AssemblyConfig) -> BreakendWriter {
        let mut w = BreakendWriter { writer: None, tumor_count: config.tumor_ids.len(), saga: config.saga_fasta_file.is_some() };
        if !config.write_types.contains(&WriteType::Breakend) {
            return w;
        }
        let mut header = vec![
            "Id", "PhaseGroupId", "PhaseSetId", "AssemblyId", "MateId", "AssemblyInfo", "Type", "Chromosome", "Position", "Orientation", "MateChr", "MatePos",
            "MateOrient", "Length", "InsertedBases", "Homology", "ConfidenceInterval", "InexactOffset", "Qual", "SplitFragments", "RefSplitFragments", "DiscFragments",
            "RefDiscFragments", "ForwardReads", "ReverseReads", "SequenceLength", "SegmentCount", "SegmentIndex", "SequenceIndex", "AlignedBases", "MapQual", "Score",
            "AdjAlignedBases", "AvgFragmentLength", "IncompleteFragments", "BreakendQual", "FacingBreakendIds", "AltAlignments", "InsertionType", "UniqueFragPos",
            "StrandBias", "ThreePrimeRange", "ClosestAssembly", "NonPrimaryFragments",
        ];
        if w.saga {
            header.push("SagaMatch");
            header.push("SagaInferred");
        }
        match std::fs::File::create(config.output_filename(WriteType::Breakend)) {
            Ok(f) => {
                let mut out = BufWriter::new(f);
                let _ = writeln!(out, "{}", header.join("\t"));
                w.writer = Some(out);
            }
            Err(e) => crate::error!("failed to initialise breakend writer: {}", e),
        }
        w
    }

    pub fn close(&mut self) {
        if let Some(w) = self.writer.as_mut() {
            let _ = w.flush();
        }
    }

    pub fn write_breakends(&mut self, assembly_alignment: AsmAlignP) {
        let Some(writer) = self.writer.as_mut() else { return };
        let assembly_info = assembly_alignment.info();
        let mut closest_assembly_map: HashMap<usize, String> = HashMap::new();
        let assemblies = assembly_alignment.assemblies();
        let first_assembly = assemblies[0];
        // a null pointer exception in Java if the first assembly has no phase group and there are breakends
        let breakends = assembly_alignment.breakends();
        for breakend in breakends {
            let phase_group = first_assembly.phase_group().expect("phase group");
            let mut sj: Vec<String> = Vec::with_capacity(48);
            sj.push(breakend.id().to_string());
            sj.push(phase_group.id().to_string());
            sj.push(assembly_alignment.phase_set().map(|p| p.id()).unwrap_or_else(|| first_assembly.non_phase_set_id()).to_string());
            sj.push(assembly_alignment.id().to_string());
            sj.push(breakend.other_breakend().map(|o| o.id().to_string()).unwrap_or_default());
            sj.push(assembly_info.clone());
            sj.push(breakend.sv_type().name().to_string());
            sj.push(breakend.chromosome.clone());
            sj.push(breakend.position.to_string());
            sj.push(breakend.orient.to_string());
            if let Some(other) = breakend.other_breakend() {
                sj.push(other.chromosome.clone());
                sj.push(other.position.to_string());
                sj.push(other.orient.to_string());
                sj.push(breakend.sv_length().to_string());
            } else {
                sj.extend(["", "", "", "0"].map(String::from));
            }
            sj.push(breakend.inserted_bases.clone());
            if breakend.homology.exists() {
                sj.push(breakend.homology.homology.clone());
                sj.push(format!("{},{}", breakend.homology.exact_start, breakend.homology.exact_end));
                sj.push(format!("{},{}", breakend.homology.inexact_start, breakend.homology.inexact_end));
            } else {
                sj.extend(["", "0,0", "0,0"].map(String::from));
            }
            sj.push(breakend.calc_sv_qual().to_string());
            let (mut split_frags, mut ref_split_frags, mut disc_frags, mut ref_disc_frags, mut forward_reads, mut reverse_reads) = (0, 0, 0, 0, 0, 0);
            for (i, support) in breakend.sample_support().iter().enumerate() {
                split_frags += support.split_fragments;
                disc_frags += support.discordant_fragments;
                if i >= self.tumor_count {
                    ref_split_frags += support.split_fragments;
                    ref_disc_frags += support.discordant_fragments;
                }
                forward_reads += support.forward_reads;
                reverse_reads += support.reverse_reads;
            }
            for v in [split_frags, ref_split_frags, disc_frags, ref_disc_frags, forward_reads, reverse_reads, assembly_alignment.full_sequence_length()] {
                sj.push(v.to_string());
            }
            let segment = &breakend.segments()[0];
            sj.push(breakend.segments().len().to_string());
            sj.push(segment.index.to_string());
            sj.push(segment.sequence_index.to_string());
            sj.push(segment.alignment.aligned_bases().to_string());
            sj.push(segment.alignment.map_qual().to_string());
            sj.push(segment.alignment.score().to_string());
            sj.push(segment.alignment.adjusted_alignment().to_string());
            sj.push(breakend.average_fragment_length().to_string());
            sj.push(breakend.incomplete_fragment_count().to_string());
            sj.push(breakend.calc_qual().to_string());
            sj.push(breakend.facing_breakends().iter().map(|x| x.id().to_string()).collect::<Vec<_>>().join(";"));
            sj.push(AlternativeAlignment::to_vcf_tag(breakend.alternative_alignments()));
            sj.push(get_insertion_type(*breakend, &assemblies, breakends).name().to_string());
            sj.push(breakend.unique_fragment_position_count().to_string());
            sj.push(breakend.sample_support().iter().map(|x| java_fixed(x.strand_bias(), 2)).collect::<Vec<_>>().join(";"));
            let read_orientation_range = breakend.read_orientation_range();
            sj.push(format!("{};{}", read_orientation_range[0], read_orientation_range[1]));
            sj.push(get_closest_assembly(*breakend, &assemblies, &mut closest_assembly_map, true));
            sj.push(breakend.non_primary_assembly_fragment_count().to_string());
            if self.saga {
                sj.push(assembly_alignment.saga_match().map(|m| m.variant_string()).unwrap_or_default());
                sj.push(breakend.is_saga_inferred().to_string());
            }
            if let Err(e) = writeln!(writer, "{}", sj.join("\t")) {
                crate::error!("failed to write breakends: {}", e);
            }
        }
    }
}

fn get_insertion_type(breakend: BreakendP, assemblies: &[AsmP], breakends: &[BreakendP]) -> InsertionType {
    if !assemblies.iter().any(|x| x.has_line_sequence()) {
        return InsertionType::None;
    }
    if is_mobile_line_element(breakend.orient, &breakend.inserted_bases) {
        return InsertionType::Line;
    }
    for other_breakend in breakends {
        if *other_breakend == breakend {
            continue;
        }
        if within_line_proximity(breakend.position, other_breakend.position, breakend.orient, other_breakend.orient)
            && is_mobile_line_element(other_breakend.orient, &other_breakend.inserted_bases)
        {
            return InsertionType::Near;
        }
    }
    InsertionType::None
}

fn get_closest_assembly(breakend: BreakendP, assemblies: &[AsmP], closest_assembly_map: &mut HashMap<usize, String>, check_mate: bool) -> String {
    if let Some(s) = closest_assembly_map.get(&breakend.addr()) {
        return s.clone();
    }
    for assembly in assemblies {
        let junction = assembly.junction();
        if junction.chromosome == breakend.chromosome && junction.orient == breakend.orient && (junction.position - breakend.position).abs() < 100 {
            let junction_coords = junction.coords_typed();
            if !breakend.is_single() {
                closest_assembly_map.insert(breakend.addr(), junction_coords.clone());
            }
            return junction_coords;
        }
    }
    let Some(other_breakend) = breakend.other_breakend() else { return String::new() };
    closest_assembly_map.insert(breakend.addr(), String::new());
    if !check_mate {
        return String::new();
    }
    let other_breakend_coords = match closest_assembly_map.get(&other_breakend.addr()) {
        Some(s) => s.clone(),
        None => {
            let s = get_closest_assembly(other_breakend, assemblies, closest_assembly_map, false);
            closest_assembly_map.insert(other_breakend.addr(), s.clone());
            s
        }
    };
    if !other_breakend_coords.is_empty() { format!("{other_breakend_coords}{OTHER_BREAKEND_COORDS}") } else { String::new() }
}
