# SIGS in Rust: record of the port

This report covers how `sigs-rs` was checked against Java SIGS 1.3, the faults found in the Java code, and the
timings. Java SIGS 1.3 is this repository's `sigs/target/sigs-1.3-jar-with-dependencies.jar`, built from sources that follow upstream commit ca04e6c39d
of 2026-11-20. Options and outputs are in [README.md](README.md).

## How it was checked

The port ran with `-java_compat` against the jar (Java 23, `-Xmx8g`; the program has no thread option) on the
same inputs, one run at a time on a shared Mac. "Identical" means every file Java wrote equals the port's byte
for byte, the same files exist, and the exit codes agree. Where Java fails, the port fails on the same input
and leaves the same files written. Times are wall times (a few tens of milliseconds are noise at a load
average of about 4), and memory is the peak resident set (`/usr/bin/time -l`).

Inputs:
- **GIAB HG008** (public), from the oncoanalyser run of 2026-11-26, with the pipeline's command line and
  resources.
- **A private whole-genome pair** from an all-Java oncoanalyser run, as renamed copies (TUMOR / NORMAL).
- **Synthetic inputs** made for the port, with no sample data in them.

**Against the pipeline's own outputs** (`sigs:1.3-beta.2`): for HG008 and the private pair, the jar's outputs equal the
pipeline's, file for file.

## What is ported

All of `SignatureFitter`, the program oncoanalyser runs, with every option.

| Java | Rust |
|---|---|
| SignatureFitter, SigsConfig, CommonUtils | `fitter.rs`, `config.rs`, `main.rs` |
| VariantFilters, SigVariant, SigSnvLoader; hmf-common SmallVariantFactory (the parts SIGS reads), SnvSigUtils, HumanChromosome.contains, Nucleotides.swapDnaBase | `snv_loader.rs`, `config.rs` |
| hmf-common LeastSquaresFit | `lsq.rs` |
| hmf-common Matrix, MatrixFile | `matrix.rs` |
| hmf-common SigUtils, SigResiduals, DataUtils.round/sizeToStr, SignatureAllocationFile | `fitter.rs` |
| hmf-common PositionFrequencies; PositionFreqBuilder.writeSampleCounts | `position_freq.rs` |
| hmf-common resource `sigs/snv_cosmic_signatures.csv` (the default signatures) | `resources/snv_cosmic_signatures.csv` (identical to `data/resources/misc/sigs/`) |

Not ported: `PositionFreqBuilder.main`, a second entry point that builds cohort position-frequency matrices
(oncoanalyser does not run it), and `SigIndelLoader` and `SigMnvLoader`, which nothing in the jar calls.

## Java behaviour reproduced (not faults)

- **The fit:** LeastSquaresFit (Lawson and Hanson's NNLS as translated by Kaminsky) step for step, including
  its reuse of the loop variable `jj` across the back substitutions; no fused multiply-adds, so the
  contributions are Java's bits.
- **Rounding:** `DataUtils.round` is `Math.round(v / 10^-n) * 10^-n` with Java's `Math.pow(10, -n)`. For n = 5
  that is `0x3ee4f8b588e368f0`, one unit in the last place below the literal `1e-5` (checked with a Java
  harness); the port takes StrictMath's pow from hmf-common-rs, which gives the same bits for 10^-3 and
  10^-5.
- **The variant factory:** a record is used when it is not filtered, on 1-22/X/Y (`chr` optional), without
  an NT other than `ref`, with allelic depths of at least two values and more than zero reads (DP or the AD
  sum, whichever is larger). An invalid `SOM_LH` ends the run for any such record, SNV or not, as in Java.
- **Records Java cannot count end the run:** a context with `N` is skipped, but a missing or lower-case
  context gives a bucket name that does not exist (Java logs it and exits 1), or, for an A or G reference,
  a StringIndexOutOfBoundsException. SAGE writes upper-case contexts, so neither occurs in pipeline data.
- **A signatures file of other than 96 rows:** the counts are not copied into the fit (Java's copyVector
  skips arrays of another length), and the residuals end the run with an ArrayIndexOutOfBoundsException.
- **Multi-sample output** uses Java's `%g`, added to hmf-common-rs as `fmt::java_general` (trailing zeros kept).
- **Registered defaults that are never used:** `-subclonal_max` and `-ploidy_max` are registered with a default
  of 1.0 (and described as "Subclonal minimum" and "Ploidy minimum"), but the filters apply only to values
  given on the command line, as in Java.

## Cases and results

All against Java with `-java_compat`; `t0` is the private pair. "files" counts the files compared.

| Case | Java time, memory | Rust time, memory | Exit | Result |
|---|---|---|---|---|
| hg008 (pipeline command) | 0.72 s, 332 MB | 0.11 s, 4.3 MB | 0/0 | identical (2 files) |
| t0 (pipeline command) | 0.54 s, 248 MB | 0.04 s, 3.1 MB | 0/0 | identical (2) |
| hg008, default signatures | 0.71 s, 356 MB | 0.11 s, 4.3 MB | 0/0 | identical (2) |
| hg008, gzipped signatures file | 0.70 s, 354 MB | 0.11 s, 4.3 MB | 0/0 | identical (2) |
| hg008 `-tumor_only` | 0.71 s, 359 MB | 0.12 s, 4.7 MB | 0/0 | identical (2) |
| t0 `-tumor_only` | 0.54 s, 251 MB | 0.04 s, 3.0 MB | 0/0 | identical (2) |
| hg008 `-ploidy_min 0.5 -ploidy_max 1.5 -subclonal_max 0.5` | 0.73 s, 353 MB | 0.12 s, 4.6 MB | 0/0 | identical (2) |
| hg008 `-subclonal_min 0.5` | 0.71 s, 360 MB | 0.12 s, 4.7 MB | 0/0 | identical (2) |
| hg008 `-min_alloc 100 -min_alloc_perc 0.05` | 0.70 s, 328 MB | 0.11 s, 4.1 MB | 0/0 | identical (2) |
| hg008 `-output_id test1` | 0.70 s, 355 MB | 0.11 s, 4.3 MB | 0/0 | identical (2) |
| hg008 `-position_bucket_size 1000000 -ref_genome_version 38` (SI1) | 0.71 s, 359 MB | 0.12 s, 4.3 MB | 0/0 | identical (3) |
| the same with `-write_pos_freq_coords` | 0.72 s, 359 MB | 0.12 s, 4.4 MB | 0/0 | identical (3) |
| bucket size 300,000,000 with coordinates (SI2) | 0.70 s, 359 MB | 0.11 s, 4.2 MB | 0/0 | identical (3) |
| bucket size 500,000, default version 37 | 0.72 s, 367 MB | 0.12 s, 4.4 MB | 0/0 | identical (3) |
| two samples by `-sample_id_file` (tumor and normal of HG008) | 1.01 s, 695 MB | 0.22 s, 5.8 MB | 0/0 | identical (1, `%g`) |
| no `-output_dir` (SI5) | 0.75 s, 364 MB | 0.12 s, 4.5 MB | 0/0 | identical (2) |
| a sample not in the VCF (SI3) | 0.60 s, 249 MB | 0.12 s, 3.1 MB | 0/0 | identical (2) |
| no sample (SI7) | 0.35 s, 95 MB | 0.00 s, 1.8 MB | 0/0 | identical (none written) |
| no VCF or counts file (SI7) | 0.35 s, 96 MB | 0.00 s, 1.8 MB | 0/0 | identical (none) |
| synthetic, 4,815 records (below) | 0.46 s, 195 MB | 0.01 s, 3.1 MB | 0/0 | identical (2) |
| the same as a plain VCF, `-log_debug` | 0.50 s, 196 MB | 0.02 s, 3.1 MB | 0/0 | identical (2) |
| synthetic, `-tumor_only -ploidy_min 1 -subclonal_max 0.8` | 0.48 s, 193 MB | 0.01 s, 3.2 MB | 0/0 | identical (2) |
| synthetic, the normal sample | 0.47 s, 194 MB | 0.02 s, 3.2 MB | 0/0 | identical (2) |
| a lower-case context | 0.34 s, 102 MB | 0.00 s, 2.2 MB | 1/1 | same failure |
| header only | 0.36 s, 104 MB | 0.00 s, 2.1 MB | 0/0 | identical (2) |
| no signature with counts in the sample's buckets (SI6) | 0.36 s, 104 MB | 0.00 s, 2.2 MB | 0/0 | identical (2) |
| a signatures file of 95 rows | 0.70 s, 351 MB | 0.11 s, 4.2 MB | 1/1 | same failure, counts file identical |
| no context, C reference | 0.35 s, 103 MB | 0.01 s, 2.2 MB | 1/1 | same failure |
| no context, A reference | 0.37 s, 104 MB | 0.01 s, 2.2 MB | 1/1 | same failure |
| an unknown SOM_LH | 0.35 s, 103 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| `-sample_counts_file` with `-sample` (SI4) | 0.32 s, 96 MB | 0.00 s, 2.0 MB | 0/0 | identical (1) |
| `-sample_counts_file` with `-sample_id_file` (SI4) | 0.33 s, 97 MB | 0.00 s, 2.0 MB | 0/0 | identical (1) |

The **synthetic VCF** has 4,800 seeded random SNVs over chr1-22, X and Y (random contexts, depths, SUBCL,
PURPLE_VCN, some SOM_LH, HOTSPOT, GND_FREQ and PON-filtered), plus edge records:
- NT `het` and `ref`, a context with N;
- a multi-allelic SNV, an MNV, an insertion;
- no AD, no reads, DP below the AD sum, lower-case alleles;
- an unplaced contig, chrM, a contig without `chr`;
- FILTER `.` and `PASS;PON`.

It is not sorted, so it also tests reading a VCF without an index.

**The corrected (default) mode against Java:** the same files in 21 of the 32 cases. The others:
- the four position-frequency cases: the counts file (SI1, SI2);
- no `-output_dir`: the files' names (SI5);
- a sample not in the VCF, no sample, no VCF: exit status 1 (SI3, SI7);
- the zero fit: the allocation file (SI6);
- the two counts-file cases (SI4).

On the pipeline inputs (HG008 and the private pair) the corrected mode gives Java's files.

## Findings

Each is corrected by default behind `compat::java("SIn")`; `-java_compat` gives Java's behaviour.

- **SI1 (SigSnvLoader / PositionFrequencies): `-position_bucket_size` always writes zeros.**
  - **What Java does:** the loader builds its PositionFrequencies as `new PositionFrequencies(V37, size)`,
    whose constructor sets `buildMap = false`. `addPosition` then never fills the per-bucket map that
    `PositionFreqBuilder.writeSampleCounts` reads.
  - **Effect:** `<sample>.sig.pos_freq_counts.csv` has a 0 for every bucket, whatever the VCF.
  - **Evidence:** HG008 with a 1 Mb bucket size: 3,102 buckets, none non-zero, for 9,227 SNVs counted.
  - **Corrected:** the map is kept. The writer matches its chromosome names to those of
    `-ref_genome_version` with or without "chr", so the default (37) works on GRCh38 names too, and adds up
    the counts. The bucket totals then equal the SNVs counted (9,227 in each HG008 case).
  - **Off by default:** oncoanalyser does not pass the option.
- **SI2 (PositionFrequencies.initialisePositionCache): chromosomes after a one-bucket first chromosome all
  start at bucket 0.**
  - **What Java does:** the start of the next chromosome is `last + 1` only when `last > 0`. A bucket size of at
    least chr1's length makes chr1 one bucket (last = 0), and every chromosome after it starts at 0 too.
  - **Effect:** all counts land in bucket 0, and in the coordinates file every bucket is named chrY. Its
    positions are `index * size`, which overflows Java's `int` ("chrY,-1989934592").
  - **Corrected:** each chromosome starts after the previous one; the usual bucket sizes are unaffected.
- **SI3 (SmallVariantFactory.createVariant): a sample that is not in the VCF gives empty results.** Java
  finds no genotype, so it makes no variant. It writes zero counts and an allocation file with a header
  only, and exits 0. Corrected: an error naming the sample, exit 1.
- **SI4 (CommonUtils.loadSampleMatrixCounts): `-sample_counts_file` cannot work.**
  - **What Java does:** it passes the configured sample list as the list to which MatrixFile adds the
    file's column names. So the list grows by the header (S1 becomes S1, S1, S2), and BucketName is looked
    for at its index in that longer list.
  - **Effect:** every row then has one field too few, so it is "invalid" and stays 0. Java fits three
    samples of zero counts and writes `sig_snv_fit.csv` with a header only. With other headers the shifted
    index makes Java parse a bucket name as a number and fail.
  - **Corrected:** the file's columns are read into a list of their own, and each configured sample takes the
    column of its name (an error if it has none).
  - **Intent:** that the configured samples select columns is my reading of the intent; Java never got that
    far.
- **SI5 (FileWriterUtils.parseOutputDir with no `-output_dir`): file names start with "null".** Java joins the
  null directory into the names. It writes `nullHG008-T.sig.allocation.tsv` and
  `nullHG008-T.sig.snv_counts.csv` in the working directory. Corrected: the working directory, without the
  prefix.
- **SI6 (SignatureFitter.processSampleResults): a fit of zero is written as zero allocations.**
  - **What Java does:** where no signature has counts in the sample's buckets, the fit is all zeros. Its
    scaling by `total / 0` makes every allocation NaN. NaN passes both minimums, rounds to 0, and is written.
  - **Effect:** every signature at 0.000, and MISALLOC at 0.000, though nothing is allocated.
  - **Corrected:** no scaling. No signature is written, and MISALLOC holds all the counts (50.000, 1.00000 in
    the case).
  - **Cannot happen with COSMIC:** every bucket of the COSMIC set has some weight.
- **SI7 (SignatureFitter.initialise): configuration errors exit with status 0.** "no sample ID(s)
  configured" and "no sample VCF or sample count file configured" are logged, and Java returns normally. A
  pipeline step would then fail later on missing files. Corrected: exit status 1.

## Speed and memory

6.5 (HG008) and 13 (private pair) times faster than Java on the pipeline inputs, in 3-4 MB against 250-330 MB.

Where the time goes (HG008, 0.11 s):
- inflating the 3.6 MB somatic VCF (17 MB of text): 0.02 s;
- reading its 17,254 lines and parsing the 10,839 unfiltered ones: 0.08 s. The filtered records and other
  contigs are skipped before parsing, using hmf-common-rs's new `VcfReader::next_line` (with every record
  parsed the run took 0.15 s);
- the fit (96 buckets by 30 signatures) and writing: under 1 ms.

Java's 0.7 s is mostly JVM start and htsjdk's decoding of every record into a full SmallVariant. SIGS has no
thread option in Java or the port.

## Open points

- GRCh37 inputs are not tested (no 37 VCF at hand). Only the position buckets depend on the version.
- Records are parsed as hmf-common-rs's `variant_context` parses them. A filtered record that htsjdk would
  reject as malformed is skipped by the port without a look, where Java fails.

## Shared code

Added to `hmf-common-rs` for this port (additions only; no existing behaviour changed):

- `fmt::java_general(d, prec)`: `String.format("%.<prec>g", d)` as Java writes it:
  - half-up to `prec` significant digits;
  - plain for decimal exponents -4 to `prec - 1`, else scientific with a two-digit exponent;
  - trailing zeros kept;
  - NaN and Infinity.

  Test `general_tests::java_percent_g`: 22 values (and their negatives) checked against Java 23's output.
  Used by SIGS's multi-sample fit file.
- `variant_context::VcfReader::next_line()`: the next data line as read, unparsed, for readers that check a
  few columns before `VcfRecord::parse` (SIGS skips filtered records this way; PEACH skips records off the
  panel).

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 9 on Google Batch ran SIGS, Virus Interpreter, PEACH and CHORD as ports, in Java mode, with every
other step from the first run's cache. SIGS's outputs (both files) are identical to the first run's Java
outputs. CUPPA and ORANGE, which ran again on them, gave the same outputs as before. SIGS took 574 ms, 2.8 MB,
against Java's 4 s, 266 MB. The rerun's tasks cost $0.02.
