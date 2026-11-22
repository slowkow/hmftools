//! ESVEE assembly step (Java package esvee.assembly, class AssemblyApplication).
#![allow(dead_code)]

pub mod alignment;
pub mod application;
pub mod assembly_deduper;
pub mod assembly_utils;
pub mod common_utils;
pub mod config;
pub mod constants;
pub mod discordant_junction_builder;
pub mod extension_seq_builder;
pub mod globals;
pub mod indel_builder;
pub mod junction_assembler;
pub mod junction_group_assembler;
pub mod junction_read_types;
pub mod line_utils;
pub mod mem_stats;
pub mod orientation;
pub mod perf;
pub mod output;
pub mod phase;
pub mod ptr;
pub mod read;
pub mod read_parse_state;
pub mod ref_base_extender;
pub mod ref_base_seq_builder;
pub mod ref_genome;
pub mod remote_region_finder;
pub mod saga;
pub mod sam;
pub mod seq_tech_utils;
pub mod sequence_builder;
pub mod sequence_compare;
pub mod sequence_diff;
pub mod types;
pub mod vis;

#[cfg(test)]
mod tests;

use anyhow::Result;

pub fn run(argv: &[String]) -> Result<()> {
    application::run(argv)
}
