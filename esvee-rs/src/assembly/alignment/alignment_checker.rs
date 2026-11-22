//! AlignmentChecker: tests an assembly against a decoy genome, and realigns the ref bases of assemblies
//! whose reads have low mapping quality.

use std::io::Write;
use std::sync::{Arc, Mutex};

use crate::assembly::alignment::bwa::{BwaMemAlignment, IBwaMemAligner};
use crate::assembly::config::AssemblyConfig;
use crate::assembly::constants::*;
use crate::assembly::sam::{cigar_elements_from_str, CigarOp};
use crate::assembly::types::junction_assembly::JunctionAssembly;

pub type TsvWriter = Arc<Mutex<std::io::BufWriter<std::fs::File>>>;

pub struct AlignmentChecker {
    assembly_map_qual_threshold: i32,
    decoy_aligner: Option<Arc<dyn IBwaMemAligner>>,
    aligner: Option<Arc<dyn IBwaMemAligner>>,
    writer: Option<TsvWriter>,
    /// decoy match rows not yet written: the application writes them in junction group order, which is the
    /// order of a single-threaded Java run
    pending_rows: std::cell::RefCell<Vec<String>>,
}

impl AlignmentChecker {
    /// `aligner` is the reference genome aligner with default parameters (used when the map-qual threshold
    /// is above zero) and `decoy_aligner` the decoy genome aligner, if configured.
    pub fn new(config: &AssemblyConfig, aligner: Option<Arc<dyn IBwaMemAligner>>, decoy_aligner: Option<Arc<dyn IBwaMemAligner>>, writer: Option<TsvWriter>) -> AlignmentChecker {
        AlignmentChecker {
            assembly_map_qual_threshold: config.assembly_map_qual_threshold,
            decoy_aligner,
            aligner: if config.assembly_map_qual_threshold > 0 { aligner } else { None },
            writer,
            pending_rows: std::cell::RefCell::new(Vec::new()),
        }
    }

    /// The decoy match rows formed since the last call.
    pub fn take_pending_rows(&self) -> Vec<String> {
        std::mem::take(&mut *self.pending_rows.borrow_mut())
    }

    pub fn matches_decoy(&self, assembly: &JunctionAssembly) -> bool {
        let Some(decoy_aligner) = &self.decoy_aligner else { return false };
        let full_sequence = assembly.bases();
        let alignment_results = decoy_aligner.align_sequence(full_sequence);
        if alignment_results.is_empty() {
            return false;
        }
        let top_alignment = &alignment_results[0];
        if top_alignment.aligner_score == 0 || top_alignment.cigar.is_empty() {
            return false;
        }
        let next_alignment = alignment_results.get(1);
        let aligned_bases: i32 = cigar_elements_from_str(&top_alignment.cigar).iter().filter(|x| x.op == CigarOp::M).map(|x| x.len).sum();
        if self.writer.is_some() {
            self.pending_rows.borrow_mut().push(match_info_row(assembly, top_alignment, aligned_bases, next_alignment));
        }
        let sequence_length = full_sequence.len() as i32;
        let unaligned_bases = sequence_length - aligned_bases;
        let score_factor = top_alignment.aligner_score as f64 / sequence_length as f64;
        unaligned_bases <= DECOY_MAX_MISMATCHES && score_factor >= DECOY_MIN_SCORE_FACTOR
    }

    pub fn fails_mappability(&self, assembly: &JunctionAssembly) -> bool {
        let Some(aligner) = &self.aligner else { return false };
        if self.assembly_map_qual_threshold == 0 || average_assembly_support_map_qual(assembly) >= self.assembly_map_qual_threshold as f64 {
            return false;
        }
        let ref_base_sequence = assembly.form_ref_base_sequence();
        let alignment_results = aligner.align_sequence(ref_base_sequence.as_bytes());
        if alignment_results.is_empty() {
            return true;
        }
        for alignment in &alignment_results {
            if alignment.map_qual >= self.assembly_map_qual_threshold {
                return false;
            }
            if alignment.xa_tag.is_some() {
                return false;
            }
        }
        true
    }
}

fn average_assembly_support_map_qual(assembly: &JunctionAssembly) -> f64 {
    let total: i32 = assembly.support().iter().map(|x| x.map_qual()).sum();
    if assembly.support_count() > 0 { total as f64 / assembly.support_count() as f64 } else { 0.0 }
}

pub fn initialise_writer(config: &AssemblyConfig) -> Option<TsvWriter> {
    if !config.write_types.contains(&crate::assembly::config::WriteType::DecoyMatches) {
        return None;
    }
    config.output_dir.as_ref()?;
    let filename = config.output_filename(crate::assembly::config::WriteType::DecoyMatches);
    match std::fs::File::create(&filename) {
        Ok(f) => {
            let mut w = std::io::BufWriter::new(f);
            let header = ["AssemblyChromosome", "AssemblyPosition", "AssemblyOrientation", "SequenceLength", "AlignCigar", "AlignMatchedBases", "AlignMismatches", "AlignLocation", "AlignScore", "AlignNextScore"];
            let _ = writeln!(w, "{}", header.join("\t"));
            Some(Arc::new(Mutex::new(w)))
        }
        Err(e) => {
            crate::error!("failed to initialise decoy match writer: {}", e);
            None
        }
    }
}

fn match_info_row(assembly: &JunctionAssembly, top_alignment: &BwaMemAlignment, aligned_bases: i32, next_alignment: Option<&BwaMemAlignment>) -> String {
    let junction = assembly.junction();
    [
        junction.chromosome.to_string(),
        junction.position.to_string(),
        junction.orient.to_string(),
        (assembly.extension_length() + assembly.ref_base_length()).to_string(),
        top_alignment.cigar.clone(),
        aligned_bases.to_string(),
        top_alignment.n_mismatches.to_string(),
        top_alignment.xa_tag.clone().unwrap_or_else(|| "null".to_string()),
        top_alignment.aligner_score.to_string(),
        next_alignment.map(|a| a.aligner_score).unwrap_or(-1).to_string(),
    ]
    .join("\t")
}

/// Writes decoy match rows to the decoy match TSV.
pub fn write_match_rows(writer: Option<&TsvWriter>, rows: &[String]) {
    let Some(writer) = writer else { return };
    let mut w = writer.lock().unwrap();
    for row in rows {
        let _ = writeln!(w, "{row}");
    }
}
