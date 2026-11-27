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
//! - `ref_coords`: GRCh37 and GRCh38 chromosome lengths (RefGenomeCoordinates).
//! - `codon`: Codons and Nucleotides (amino acids from DNA, complements); `blosum`: BlosumMapping.
//! - `bwa` (with the `bwa` feature): BWA-MEM alignment through the vendored bwa, as gatk-bwamem-jni.
//! - `csv`: commons-csv's CSVPrinter as the Kotlin tools use it (tab-delimited, MINIMAL quoting).
//! - `kotlin_random`: kotlin.random.Random(seed), the XorWow generator.
//! - `sam_header`: SAM header lines as htsjdk writes them (`@HD VN:1.6` whatever the input's version).
//! - `gz`: fast gzip / BGZF text reading (zlib-rs) with a line reader.
//! - `delim_file`, `redux_files`, `cobalt_files`, `esvee_files`: DelimFileReader, the plain or gzip
//!   writer of FileWriterUtils.createBufferedWriter and DelimFileWriter's rows; REDUX's BQR, duplicate
//!   frequency and jitter count files, COBALT's GC median file and ESVEE's discordant statistics (qsee).
//! - `vcf`: VCF headers, records and the writer as htsjdk 2.23 has them, with the BGZF stream and the
//!   tabix index written as Java writes them.
//! - `variant_context`: VCF records with htsjdk's typed attributes and lazy genotypes (VariantContext),
//!   for tools that read a VCF, change fields and write it again (PURPLE); written through `vcf`'s BGZF
//!   stream and tabix index.
//! - `commons_math`, `java_strictmath`, `kernel_estimator`, `config_builder`: commons-math3
//!   distributions, StrictMath, KernelEstimator and ConfigBuilder as the Java tools use them.
//! - `neo`, `hla`, `rna`: the neoepitope files of LINX, NEO and ISOFOX, LILAC's allele file, and ISOFOX's
//!   transcript expression and the cohort expression matrix (NEO).

pub mod args;
pub mod blosum;
#[cfg(feature = "bwa")]
pub mod bwa;
pub mod codon;
pub mod cobalt_files;
pub mod commons_math;
pub mod config_builder;
pub mod compat;
pub mod csv;
pub mod delim_file;
pub mod driver_panel;
pub mod ensembl;
pub mod esvee_files;
pub mod faidx;
pub mod fmt;
pub mod gz;
pub mod hla;
pub mod java_fastmath;
pub mod java_hash;
pub mod java_id_hash;
pub mod java_math;
pub mod java_sort;
pub mod java_strictmath;
pub mod java_text;
pub mod kernel_estimator;
pub mod known_fusion;
pub mod kotlin_random;
pub mod log;
pub mod metrics;
pub mod neo;
pub mod ref_coords;
pub mod refgenome;
pub mod redux_files;
pub mod sam_header;
pub mod regions;
pub mod rna;
pub mod poisson;
pub mod purple_files;
pub mod segmentation;
pub mod sigs;
pub mod sv;
pub mod transcript_utils;
pub mod variant_context;
pub mod vcf;
