# teal-rs

A port of hmftools TEAL 1.4 (Kotlin) to Rust: telomeric content and length, and telomeric
rearrangements (break ends) from whole-genome BAMs. It takes the same command lines as the Java tool and
writes the same files; with `-java_compat` the outputs equal Java's, byte for byte after decompression.
Without it, the Java faults found in the port (TE1 to TE11 in [REPORT.md](REPORT.md)) are corrected.

## Usage

```
cargo build --release          # target/release/teal-rs
teal-rs [<class> | prep | pipeline | app] <options>
```

The image's `teal` command (`../rust-docker/bin/teal`) runs `teal-rs` with oncoanalyser's command lines.

| Class (short name) | What it does | Outputs |
|---|---|---|
| `com.hartwig.hmftools.teal.TealPipelineTelbamApp` (`prep`) | reads each BAM or CRAM and collects the telomeric reads with all reads of their fragments | `<id>.teal.telbam.bam` per sample |
| `com.hartwig.hmftools.teal.TealPipelineApp` (`pipeline`) | from the telbams and the COBALT, PURPLE and BamTools outputs: telbams again, telomere length per sample, break ends for a tumor/normal pair | `<id>.teal.telbam.bam`, `<id>.teal.tellength.tsv`, `<tumor>.teal.breakend.tsv.gz` |
| `com.hartwig.hmftools.teal.TealApplication` (`app`, also the default) | both, with the read depths on the command line | as above |

JVM options (`-Xmx...`, `-Xms...`, `-XX:...`, `-D...=...`) are accepted and ignored, as oncoanalyser calls
the bioconda wrapper with them. `teal-rs -version` prints `Teal version 1.4`; `-help` lists the options.

Common options (all entry points): `-reference`, `-reference_bam`, `-tumor`, `-tumor_bam` (one sample of
each, or one alone: germline-only or tumor-only), `-ref_genome` (for CRAM), `-ref_genome_version` (37 or
38, default V37), `-sequencing_type` (ILLUMINA, ULTIMA, SBX), `-threads`, `-output_dir`, `-log_debug`,
`-log_level`.

TealPipelineApp: `-cobalt_dir` (`<id>.cobalt.gc.median.tsv`), `-purple_dir` (`<tumor>.purple.purity.tsv`
and `.purple.qc`), `-reference_wgs_metrics` and `-tumor_wgs_metrics` (BamTools `bam_metric.summary.tsv`).
TealApplication: `-reference_mean_read_depth`, `-reference_gc50_read_depth`,
`-reference_duplicate_proportion`, `-tumor_mean_read_depth`, `-tumor_gc50_read_depth`,
`-tumor_duplicate_proportion`, `-tumor_purity`, `-tumor_ploidy`.

Port options: `-java_compat` (or `TEAL_JAVA_COMPAT=1`) reproduces every Java fault;
`-java_compat_except TE2,TE4` (or `TEAL_JAVA_COMPAT_EXCEPT`) all but those listed.

## How it works

The prep reads every record of the BAM. Java does so in partitions (the unplaced unmapped reads first,
then 10 Mb regions in header order) on `threads - 1` reader threads; with one reader thread (`-threads` 1
or 2) its result is that of reading the records in that order, which the port does for every thread count.
BAMs are read with the port's own BGZF reader (`src/bgzf.rs`): batches of blocks inflated by every thread
(libdeflate), records parsed in place; only records that are telomeric or belong to an open read group
become htslib records. The telomeric test works on the packed 4-bit sequence (SIMD on aarch64 and x86-64).
CRAM goes through htslib. The telbam is sorted as htsjdk's writer sorts it (SAMRecordCoordinateComparator,
stable) and written with the input header.

The break end analysis aligns soft clips and mates to a telomere template (hmf-common's Smith-Waterman
LocalSequenceAligner, its tie rules kept); results are cached per read, and only read groups with a read
within 500 bases of a break end are looked at.

## Tests

```
CARGO_BUILD_JOBS=2 cargo test --release
```

runs the ported Kotlin tests (TealUtilsTest, ReadGroupTest, TelomereMatcherTest, TelomereLengthCalcTest,
BreakEndSupportCounterTest, TelomericSplitReadAnalyserTest) and the port's own. The comparisons with Java
are in [REPORT.md](REPORT.md). `scripts/make_synthetic.sh <out dir> [seed] [unknown]` writes the synthetic
BAMs and pipeline files used there (needs samtools); `scripts/bam_raw_equal.py` compares two BAMs record by
record as uncompressed bytes.
