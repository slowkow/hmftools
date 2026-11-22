//! AlignmentWriter: the phased assembly TSV and the alignment TSV.

use std::io::{BufWriter, Write};

use crate::assembly::alignment::align_data::AlignP;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::common::fmt::java_fixed;

pub struct AlignmentWriter {
    saga_enabled: bool,
    writer: Option<BufWriter<std::fs::File>>,
    detailed_writer: Option<BufWriter<std::fs::File>>,
}

fn create(config: &AssemblyConfig, write_type: WriteType, header: &[&str], what: &str) -> Option<BufWriter<std::fs::File>> {
    if !config.write_types.contains(&write_type) || config.output_dir.is_none() {
        return None;
    }
    match std::fs::File::create(config.output_filename(write_type)) {
        Ok(f) => {
            let mut out = BufWriter::new(f);
            let _ = writeln!(out, "{}", header.join("\t"));
            Some(out)
        }
        Err(e) => {
            crate::error!("failed to initialise {}: {}", what, e);
            None
        }
    }
}

impl AlignmentWriter {
    pub fn new(config: &AssemblyConfig) -> AlignmentWriter {
        let saga_enabled = config.saga_fasta_file.is_some();
        let mut header = vec!["AssemblyIds", "AssemblyInfo", "Merged", "SequenceLength", "AssemblyCigar", "FullSequence"];
        if saga_enabled {
            header.push("SagaMatch");
            header.push("SagaMatchCigar");
        }
        let writer = create(config, WriteType::PhasedAssembly, &header, "alignment writer");
        let detailed_header = [
            "AssemblyIds", "AssemblyInfo", "RefInfo", "RawSeqCoords", "AdjSeqCoords", "MapQual", "Cigar", "Orientation", "AlignedBases", "Score", "Flags", "NMatches",
            "XaTag", "MdTag", "CalcAlignLength", "ModMapQual", "DroppedOnRequery", "LinkedAltAlignment",
        ];
        let detailed_writer = create(config, WriteType::Alignment, &detailed_header, "detailed alignment data writer");
        AlignmentWriter { saga_enabled, writer, detailed_writer }
    }

    pub fn close(&mut self) {
        for w in [self.writer.as_mut(), self.detailed_writer.as_mut()].into_iter().flatten() {
            let _ = w.flush();
        }
    }

    pub fn write_phased_assembly(&mut self, assembly_alignment: AsmAlignP) {
        let Some(writer) = self.writer.as_mut() else { return };
        let mut sj: Vec<String> = vec![
            assembly_alignment.assembly_ids(),
            assembly_alignment.info(),
            assembly_alignment.is_merged().to_string(),
            assembly_alignment.full_sequence_length().to_string(),
            assembly_alignment.assembly_cigar().to_string(),
            assembly_alignment.full_sequence().unwrap_or("null").to_string(),
        ];
        if self.saga_enabled {
            let saga_match = assembly_alignment.saga_match();
            sj.push(saga_match.map(|m| m.variant_string()).unwrap_or_default());
            sj.push(saga_match.map(|m| m.cigar_string()).unwrap_or_default());
        }
        if let Err(e) = writeln!(writer, "{}", sj.join("\t")) {
            crate::error!("failed to write alignment data: {}", e);
        }
    }

    pub fn write_alignment_details(&mut self, assembly_alignment: AsmAlignP, alignments: &[AlignP]) {
        let Some(writer) = self.detailed_writer.as_mut() else { return };
        if alignments.is_empty() {
            return;
        }
        let assembly_str = format!("{}\t{}", assembly_alignment.assembly_ids(), assembly_alignment.info());
        for alignment in alignments {
            let sj: Vec<String> = vec![
                assembly_str.clone(),
                alignment.ref_location_str(),
                format!("{}-{}", alignment.raw_sequence_start(), alignment.raw_sequence_end()),
                format!("{}-{}", alignment.sequence_start(), alignment.sequence_end()),
                alignment.map_qual().to_string(),
                alignment.cigar().to_string(),
                alignment.orientation().to_string(),
                alignment.aligned_bases().to_string(),
                alignment.score().to_string(),
                alignment.flags().to_string(),
                alignment.n_matches().to_string(),
                alignment.xa_tag().to_string(),
                alignment.md_tag().to_string(),
                alignment.adjusted_alignment().to_string(),
                java_fixed(alignment.modified_map_qual(), 0),
                alignment.dropped_on_requery().to_string(),
                if alignment.has_selected_alt_alignment() { alignment.selected_alt_alignment().map(|a| a.vcf_string()).unwrap_or_default() } else { String::new() },
            ];
            if let Err(e) = writeln!(writer, "{}", sj.join("\t")) {
                crate::error!("failed to write alignment details: {}", e);
            }
        }
    }
}
