//! The assembly visualiser (Java package esvee.assembly.vis): one HTML file per assembly alignment with
//! paired breakends, written to `<output_dir>/vis/` when `-run_assembly_vis` is given with specific junctions.
//!
//! Java changes objects while it draws: the cigar list of an alignment and the cigar list and base
//! qualities of a read are reversed in place, and stay that way for whatever is drawn next. The outputs of
//! the step are written before the visualiser runs, so only later drawings see this. [`VisState`] keeps
//! those lists apart from the assembly objects and carries them from one drawing to the next.

pub mod assembly_visualiser;
pub mod read_view_model;

use std::collections::HashMap;

use crate::assembly::alignment::align_data::AlignP;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::config::AssemblyConfig;
use crate::assembly::read::read::ReadP;
use crate::assembly::sam::{CigarElement, CigarOp};
use crate::common::vis::svg::VisCigarOp;

// AssemblyVisConstants
pub const READ_HEIGHT_PX: f64 = 12.0;
pub const VIS_DIR: &str = "vis";
pub const VIEW_REGION_SIZE: i32 = 500;
pub const DISPLAY_EVERY_NTH_COORD: i32 = 10;

/// The cigar list and base qualities of a read as the visualiser has left them
pub struct ReadVisState {
    pub cigar: Vec<CigarElement>,
    pub quals: Vec<u8>,
}

/// What Java's visualiser has changed in place so far, by object
#[derive(Default)]
pub struct VisState {
    align_cigars: HashMap<usize, Vec<CigarElement>>,
    reads: HashMap<usize, ReadVisState>,
}

impl VisState {
    /// AlignData.cigarElements()
    pub fn align_cigar(&mut self, alignment: AlignP) -> &mut Vec<CigarElement> {
        self.align_cigars.entry(alignment.addr()).or_insert_with(|| alignment.cigar_elements().clone())
    }
    /// Read.cigarElements() and Read.getBaseQuality()
    pub fn read(&mut self, read: &ReadP) -> &mut ReadVisState {
        self.reads.entry(read.addr()).or_insert_with(|| ReadVisState { cigar: read.cigar_elements().to_vec(), quals: read.get_base_quality().to_vec() })
    }
}

pub fn vis_cigar(elements: &[CigarElement]) -> Vec<(i32, VisCigarOp)> {
    elements
        .iter()
        .map(|e| {
            let op = match e.op {
                CigarOp::M | CigarOp::EQ | CigarOp::X => VisCigarOp::Match,
                CigarOp::S => VisCigarOp::S,
                CigarOp::I => VisCigarOp::I,
                CigarOp::D => VisCigarOp::D,
                CigarOp::H | CigarOp::N => VisCigarOp::Skip,
                CigarOp::P => VisCigarOp::Pad,
            };
            (e.len, op)
        })
        .collect()
}

/// AssemblyApplication.runVisualiser. An exception while one assembly alignment is drawn is logged at
/// debug level (Java also prints its stack trace to stderr) and the next one is tried.
pub fn run_visualiser(config: &AssemblyConfig, assembly_alignments: &[AsmAlignP]) {
    let mut state = VisState::default();
    for assembly in assembly_alignments {
        let result = match assembly_visualiser::AssemblyVisualiser::create(config, *assembly, &mut state) {
            Ok(Some(mut visualiser)) => visualiser.write_vis_file(&mut state),
            Ok(None) => Ok(()),
            Err(e) => Err(e),
        };
        if let Err(e) = result {
            if crate::common::log::debug_enabled() {
                let all = assembly_alignments.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ");
                crate::debug!("EsveeVis failed for assemblyAlignments([{}])", all);
            }
            eprintln!("{e}");
        }
    }
}

#[cfg(test)]
mod tests;
