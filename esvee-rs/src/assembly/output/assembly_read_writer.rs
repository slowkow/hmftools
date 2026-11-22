//! AssemblyReadWriter: one row per support read of each assembly.

use std::io::{BufWriter, Write};

use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::globals::{write_candidate_reads, READ_ID_TRIMMER};
use crate::assembly::sam::{FLAG_FIRST_OF_PAIR, FLAG_MATE_REVERSE_STRAND, FLAG_MATE_UNMAPPED, FLAG_READ_REVERSE_STRAND, FLAG_READ_UNMAPPED};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::SupportType;

pub struct AssemblyReadWriter {
    writer: Option<BufWriter<std::fs::File>>,
}

impl AssemblyReadWriter {
    pub fn new(config: &AssemblyConfig) -> AssemblyReadWriter {
        let mut w = AssemblyReadWriter { writer: None };
        if !config.write_types.contains(&WriteType::AssemblyRead) {
            return w;
        }
        let header = [
            "AssemblyId", "AssemblyInfo", "ReadId", "SupportType", "IsRef", "Chromosome", "PosStart", "PosEnd", "Orient", "Cigar", "InsertSize", "MateChr",
            "MatePosStart", "MatePosEnd", "Flags", "FirstInPair", "Reversed", "Unmapped", "MateUnmapped", "MateReversed", "Supplementary", "SuppData",
            "InferredFragLength", "JunctionReadStartDistance", "AlignmentIndex", "AlignmentOrientation", "BreakendSupport", "Matches", "ExtMismatches",
            "RefMismatches", "MismatchInfo", "TrimCountStart", "TrimCountEnd", "LineTail",
        ];
        match std::fs::File::create(config.output_filename(WriteType::AssemblyRead)) {
            Ok(f) => {
                let mut out = BufWriter::new(f);
                let _ = writeln!(out, "{}", header.join("\t"));
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

    pub fn write_assembly_reads(&mut self, assembly: AsmP) {
        let Some(writer) = self.writer.as_mut() else { return };
        let assembly_info = assembly.junction().coords();
        let mut support_reads: Vec<SupP> = assembly.support().clone();
        if write_candidate_reads() {
            for x in assembly.candidate_support().iter().chain(assembly.unmapped_reads().iter()) {
                support_reads.push(SupportRead::new(x, SupportType::Discordant, 0, -1, -1));
            }
        }
        let assembly_id = assembly.id().to_string();
        let mut line = String::with_capacity(512);
        for support in &support_reads {
            line.clear();
            let mut add = |s: &str| {
                if !line.is_empty() {
                    line.push('\t');
                }
                line.push_str(s);
            };
            add(&assembly_id);
            add(&assembly_info);
            add(&READ_ID_TRIMMER.restore(support.id()));
            add(support.type_().name());
            add(&support.is_reference().to_string());
            add(support.chromosome());
            add(&support.alignment_start().to_string());
            add(&support.alignment_end().to_string());
            add(&support.orientation().as_byte().to_string());
            add(support.cigar());
            add(&support.insert_size().to_string());
            add(support.mate_chromosome().unwrap_or("null"));
            add(&support.mate_alignment_start().to_string());
            add(&support.mate_alignment_end().to_string());
            add(&support.flags().to_string());
            for flag in [FLAG_FIRST_OF_PAIR, FLAG_READ_REVERSE_STRAND, FLAG_READ_UNMAPPED, FLAG_MATE_UNMAPPED, FLAG_MATE_REVERSE_STRAND] {
                add(&support.is_flag_set(flag).to_string());
            }
            add(&support.is_supplementary().to_string());
            match support.supplementary_data() {
                Some(supp) => add(&format!("{}:{}:{}", supp.chromosome, supp.position, supp.cigar)),
                None => add(""),
            }
            add(&support.inferred_fragment_length().to_string());
            add(&support.junction_read_start_distance().to_string());
            add(&support.full_assembly_index_start().to_string());
            add(&support.full_assembly_orientation().as_byte().to_string());
            add(support.breakend_support_type().map(|t| t.name()).unwrap_or("NONE"));
            add(&support.extension_base_matches().to_string());
            add(&support.extension_base_mismatches().to_string());
            add(&support.reference_mismatches().to_string());
            add(support.mismatch_info());
            add(&support.trim_count_start().to_string());
            add(&support.trim_count_end().to_string());
            add(&support.has_line_tail().to_string());
            // an empty first column cannot occur (the assembly id), so the joined line is never empty
            if let Err(e) = writeln!(writer, "{line}") {
                crate::error!("failed to write assembly reads: {}", e);
            }
        }
    }
}
