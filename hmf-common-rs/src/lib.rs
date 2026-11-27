//! Code shared by the Rust ports of hmftools: what the Java tools take from `hmf-common`, and the
//! Java library behaviour the ports reproduce.
//!
//! - `args`: command lines in the style of the Java ConfigBuilder.
//! - `compat`: the switch between corrected and Java behaviour (`-java_compat`).
//! - `faidx`, `refgenome`: indexed FASTA access with 1-based inclusive positions.
//! - `fmt`: numbers formatted as Java's `String.format` and `Double.toString` write them.
//! - `log`: log lines in the layout of the Java tools; the `info!`, `warn!`, `error!` and `debug!` macros.
//! - `regions`: chromosomes, regions and the `-specific_regions` / `-specific_chr` options.
//! - `java_hash`, `java_id_hash`, `java_sort`, `java_math`: Java's HashMap iteration order, its
//!   identity-hash order under `-XX:hashCode=2`, TimSort and stream sums, where they reach the output.
//! - `java_fastmath`: commons-math3 FastMath.log, which COBALT applies before segmenting.
//! - `segmentation`: per-arm piecewise constant segmentation and the segments file (AMBER, COBALT), and
//!   the .pcf readers (PURPLE).
//! - `ensembl`: the Ensembl data cache (genes, transcripts, exons); `driver_panel`: the driver gene panel
//!   and its gene regions; `metrics`: the BamTools metrics files (summary, coverage, fragment lengths,
//!   flag counts, gene and exon depth), writers and readers; `java_text`: Java's number parsing and
//!   String.split as the file readers need them.

pub mod args;
pub mod compat;
pub mod driver_panel;
pub mod ensembl;
pub mod faidx;
pub mod fmt;
pub mod java_fastmath;
pub mod java_hash;
pub mod java_id_hash;
pub mod java_math;
pub mod java_sort;
pub mod java_text;
pub mod log;
pub mod metrics;
pub mod refgenome;
pub mod regions;
pub mod segmentation;
