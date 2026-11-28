# Rust ports of hmftools: plan for the rest of the suite

Goal: run the whole oncoanalyser chain with Rust programs, so that more samples can be processed per
machine and per day. Each port takes the same options and writes the same files as the Java program,
corrects the faults found in Java by default, and reproduces Java exactly under `-java_compat`.

Each port started from a known upstream commit. How the ports follow upstream since then, one reviewed
change at a time, is in [`UPSTREAM.md`](UPSTREAM.md).

Done: SAGE 5.1 (`sage-rs/`), REDUX 2.1 (`redux-rs/`), ESVEE 2.1 (`esvee-rs/`). Each has a README,
PLAN and REPORT.

## 1. Where the time goes

CPU hours per step in one Illumina whole-genome tumor/normal run of oncoanalyser 2.3.0 (all Java
tools), from its Nextflow trace. CPU hours are wall time times CPU use, summed over tasks; reruns of a
step are counted once.

| Step | CPU hours | Share | Wall time | Peak memory | Status |
|---|---|---|---|---|---|
| BWA-MEM2 alignment (224 tasks) | 688.3 | 80.7% | 48.9 h summed over tasks | 25 GB | not hmftools |
| SAGE somatic | 73.3 | 8.6% | 5.0 h | 70 GB | ported |
| REDUX, tumor and normal | 35.8 | 4.2% | 8.4 h summed | 51 GB | ported |
| fastp read splitting | 17.8 | 2.1% | 1.4 h summed | 2 GB | not hmftools |
| ESVEE | 13.7 | 1.6% | 2.2 h | 121 GB | ported |
| BamTools, tumor and normal | 8.1 | 1.0% | 0.7 h summed | 4.7 GB | ported (metrics) |
| TEAL, prep and pipeline | 5.2 | 0.6% | 0.75 h | 18.5 GB | ported |
| COBALT | 4.4 | 0.5% | 0.4 h | 9.9 GB | ported |
| AMBER | 2.9 | 0.3% | 0.5 h | 9.8 GB | ported |
| VIRUSBreakend (GRIDSS) | 1.9 | 0.2% | 0.6 h | 5.9 GB | not hmftools |
| SAGE germline | 0.5 | 0.1% | 4 min | 12.7 GB | ported |
| 14 other hmftools steps | 0.7 | 0.1% | CIDER 0.9 h (3% CPU), the rest minutes | up to 10.7 GB | CIDER and LILAC ported; 12 not |

What this means:

- Alignment and read splitting are 83% of the CPU time. No hmftools port changes that.
- The three ported programs were 74% of the hmftools CPU time and set the memory needed per sample:
  Java ESVEE peaked at 121 GB, most of a 128 GB machine (Java peaks partly reflect the heap
  oncoanalyser gives each step). The ports peak at 7.2 GB (ESVEE), 1.8 GB (SAGE) and 1.6 GB (REDUX)
  on the same kind of input, so several samples can share a machine.
- The hmftools programs not yet ported used about 21 CPU hours. Most of that is five programs that
  each read a whole BAM: BamTools, TEAL, COBALT, AMBER and CIDER.
- Caveat: this run used older versions (SAGE 4.2, REDUX 1.2.2, ESVEE 1.2). oncoanalyser 3.0.0 uses
  SAGE 5.0.2, REDUX 2.0.5 and ESVEE 2.0.1; the ports follow hmftools master (5.1, 2.1, 2.1).

## 2. Before the next port: run the three ports inside oncoanalyser

The speed-ups only reach sample throughput once oncoanalyser runs the Rust programs.

1. **Command-line fit.** oncoanalyser calls each tool through its bioconda wrapper, for example
   `sage -Xmx… <options>`, and reads the version from `sage -version`. Each port must accept and
   ignore the JVM options and print its version in the same form.
2. **Builds.** Linux x86-64 release builds and one container image per tool (or one with all ports).
3. **Pipeline config.** A Nextflow config that swaps the container for the SAGE, REDUX and ESVEE
   processes, with fewer CPUs and much less memory than the Java defaults.
4. **Versions.** Run the downstream Java tools from the same hmftools commit as the ports, or check
   that PURPLE, LINX and ORANGE as pinned by oncoanalyser 3.0.0 read the newer outputs unchanged.
5. **End-to-end check** on the whole-genome pair: the Java chain against the Rust chain under
   `-java_compat` must give identical PURPLE, LINX and ORANGE outputs. Then the corrected default:
   record how purity, ploidy, drivers and structural variants change.

## 3. Alignment and read preparation

Alignment is most of the cost. Decision (2026-11-22): the pipeline will align with **bwa-mem3** in
place of BWA-MEM2. oncoanalyser calls `bwa-mem2 mem -Y -K 100000000` piped into sambamba, so the
alignment step needs its own module or container change, and the end-to-end check of section 2
should run on bwa-mem3 alignments. A Rust rewrite of the aligner is not planned.

**fastp** (2.1% of the CPU time above) changes no reads in this pipeline. oncoanalyser runs it with
adapter trimming, quality filtering, length filtering and poly-G trimming all turned off, for two
purposes only:

- splitting each FASTQ pair into chunks of `max_fastq_records` read pairs (10 million in 2.3.0), so
  that alignment can be spread over many tasks at once; each chunk is aligned and sorted separately
  and REDUX merges the chunks again;
- with UMIs enabled, moving the UMI from the read sequence into the read name (off in this run).

Its time goes into decompressing and recompressing every read. No replacement is needed:

- oncoanalyser 3.0.0 sets `max_fastq_records = -1` by default and then does not run fastp at all
  unless UMIs are enabled; in 2.3.0, `--max_fastq_records 0` does the same.
- Splitting only pays off when the chunks run on many machines at once. On one machine, a single
  aligner process per lane with all the threads avoids the split, the recompression and the merge.

Still worth measuring: starting from the provider's alignments when they are delivered on a
compatible reference, which skips alignment entirely. The hmftools callers are tuned on BWA-MEM
output (soft clips, mapping qualities, supplementary alignments), so this needs the chain run on both
alignments of the same pair and a comparison of purity, ploidy, small variants, structural variants
and drivers.

## 4. Order of the remaining ports

Sizes are lines of Java or Kotlin in `src/main`, not counting what they use from `hmf-common`. For
scale, the three ports covered about 93,000 lines (SAGE 26,700, REDUX 18,500, ESVEE 47,500).

### Wave A: programs that read whole BAMs

These run after REDUX and take 25 to 55 minutes each. AMBER and COBALT gate PURPLE; the others gate
ORANGE.

| Order | Program | Size | Why | Notes |
|---|---|---|---|---|
| 1 | AMBER | 5,800 | gates PURPLE; 2.9 CPU h | allele counts at known heterozygous sites; reuses SAGE's read handling; segmentation is Java (no R) |
| 2 | COBALT | 6,100 | gates PURPLE; 4.4 CPU h | read depth per 1 kb window with GC normalisation; segmentation is Java |
| 3 | BamTools metrics | 3,400 of 13,100 | most CPU of the unported programs, 8.1 h | only the metrics command is used by oncoanalyser |
| 4 | TEAL | 4,400 (Kotlin) | 5.2 CPU h, 18.5 GB | telomere content and rearrangements |
| 5 | CIDER | 11,500 (Kotlin) | 55 min wall at 3% CPU: waiting on reads, not computing | aligns with BWA-MEM through the same library as ESVEE, which esvee-rs already builds |
| 6 | LILAC | 13,700 | HLA typing; seconds to minutes | small runtime, needed for completeness and for NEO |

Status on 2026-11-27 (all six ported; AMBER, COBALT and BamTools also checked inside the pipeline on
GIAB HG008, where they and every later step matched the Java run):
- **AMBER** (`amber-rs`): byte-identical to Java on the whole-genome pair in all four modes. 4-7 times
  faster at 2 threads, with a fifth of the memory or less; 2.8 min against Java's 12.2 on Google Batch,
  after reading runs of regions with one index query each. Findings AM1-AM7 and AM51-AM55.
- **COBALT** (`cobalt-rs`): identical to Java on the whole-genome tumor/normal pair. 3 times faster at
  2 threads, in 0.7 GB against 4.6 GB. Findings CO1-CO6.
- **BamTools metrics** (`bam-tools-rs`): identical to Java on 19 small and option runs and on the
  whole-genome pair. 5-8 times faster (6.5 and 7.1 times on the whole-genome BAMs), in about 100 MB
  against 2-3 GB. Findings BT1-BT20.
- **TEAL** (`teal-rs`): identical to Java on 29 cases (telbams record for record, lengths, break ends);
  the prep 8-13 times faster on BAM, TealPipelineApp 19 s against 599 s on whole-genome telbams. Java
  writes some fragments twice through a race (TE1). Findings TE1-TE11. The whole-genome prep is
  identical to Java's and 11.6 times faster (649 s against 7,561 s, 1.07 GB against 2.43 GB).
- **LILAC** (`lilac-rs`): byte-identical to Java on 22 cases; under a second against 30 s, 0.4 GB
  against 5 GB. Findings LI1-LI25.
- **CIDER** (`cider-rs`): byte-identical to Java on 16 cases with BAM input (from CRAM, cider.bam's tag
  order differs). 1.5-2.3 times faster, a sixth of the memory. Java misplaces genes of ALT and patch
  contigs (CI1). Findings CI1-CI11. Its BWA-MEM bindings are in hmf-common-rs (feature `bwa`).

The pipeline's Java images for these tools are Hartwig's beta builds (`amber:4.4-beta.1`,
`cobalt:3.1-beta.1`, `bam-tools:1.7-beta.3`), from tags with the same sources as the ports follow.

Once each matches Java on its own, AMBER, COBALT, BamTools and TEAL's read collection could share a
single pass over each BAM, or be computed by REDUX while it writes the BAM. Today each program
decompresses the whole BAM again.

### Wave B: calling, annotation and reporting

Each takes minutes, so porting them saves little time; they are ported to finish the suite and remove
Java from the pipeline.

| Program | Size | Non-Java parts |
|---|---|---|
| PAVE | 9,100 | none; loads the Ensembl data cache |
| PURPLE | 22,600 | plots drawn by R scripts and Circos (Perl) |
| LINX | 41,000 | visualiser plots with R and Circos (Perl) |
| Virus Interpreter | 600 | none |
| CHORD | 2,300 | prediction by an R random forest |
| SIGS | 1,500 | none |
| PEACH | 1,900 | none |
| NEO | 15,500 | none |
| CUPPA | 3,100 Java + 10,300 Python and R | classifier in Python (scikit-learn) |
| ORANGE | 9,400 + 2,000 data model | PDF report |
| qsee | 5,800 | plots with R |
| WISP | 6,600 | plots with R |

Status on 2026-11-27 (PAVE, PURPLE and LINX are being ported first; plots are left for later, and the
pipeline's plotting steps keep their Java images):
- **PAVE** (`pave-rs`): byte-identical to Java on 29 cases (the VCF, its BGZF bytes and tabix index):
  HG008 and a private pair, somatic and germline with the pipeline's commands, slices, a synthetic set
  of 365,000 records over every kind of impact site, 13 option cases and 8 edge cases; where Java fails,
  the port fails on the same record. About 6 times faster at 2 threads, in a twentieth of the memory
  (0.2 GB against 4-5 GB). Java reads the mappability file one base off (PA2: MAPPABILITY changes on 256
  of 17,254 HG008 somatic records) and drops a gene once a variant reaches past its end (PA1). Findings
  PA1-PA9.
- **PURPLE** (`purple-rs`): byte-identical to Java (VCFs with their BGZF bytes and tabix indexes) on
  seven cases: HG008 tumor/normal, tumor-only, germline-only and targeted, and a private pair, at 1 and 2
  threads. Two outputs need a defined Java reference, as Java varies from run to run: the kataegis ids
  (PU1, a race between threads) and baf.circos (PU2, identity-hash order). 3.1-3.5 times faster
  (8.8 in germline-only mode), with 2.5-4.5 times less memory. Germline amplifications are typed as
  deletions in Java (PU9). Findings PU1-PU9. Its 9 plots (7 from R scripts, 2 from Circos) are drawn natively
  since 2026-11-27 (below); in the pipeline ORANGE runs on them, and oncoanalyser's `rust_ports` includes PURPLE.
- **LINX** (`linx-rs`): LINX_SOMATIC and LINX_GERMLINE byte-identical to Java on HG008 and on a private
  whole-genome pair (14 and 6 files) and on 14 small and option cases; 2.6-3.5 times faster at one
  thread, in a quarter to a third of the memory (0.2-0.3 GB against 0.8-1.1 GB). Java's output order
  depends on identity hashes (LX1); two of the faults are chromosome-name checks that never match on
  GRCh38 (LX2, LX5). Findings LX1-LX5. The visualiser (Circos and R in Java) is ported too and draws natively with hmf-plot-rs:
  all 5,219 data files of 13 cases byte-identical to Java's, the plots differing in rendering only; 12 s and 0.5 GB
  for HG008's 49 plots against about 5 min and 11 GB for Java in the pipeline. Findings LV1-LV4 (a 3' UTR drawn as
  a protein domain, links drawn grey). LINXREPORT (gpgr, R) is not ported.
  Inside the pipeline on HG008 (rerun 7): identical to the port run locally, and to the Java run apart
  from LX1's ids and order; ORANGE's JSON unchanged.

Second batch, merged on 2026-11-27:
- **SIGS** (`sigs-rs`), **Virus Interpreter** (`virus-interpreter-rs`) and **PEACH** (`peach-rs`): byte-identical
  to Java on 32, 11 and 25 cases (HG008, a private pair, every option, synthetic and failure cases), each in
  a tenth of a second or less and a few MB, against 0.3-0.8 s and 85-920 MB. Findings SI1-SI7 (an option that
  always writes zeros, an int overflow, configuration errors that exit 0, ...), VI1-VI3, PE1-PE2 (a record
  anywhere in the VCF that ends the run). The 1.7.1 jar calls the metrics option `-bam_metrics_tumor_dir`,
  where oncoanalyser passes `-tumor_metrics_dir`; the port takes both.

- **CHORD** (`chord-rs`): byte-identical to Java on 53 cases, the R script's randomForest prediction
  included, which the port makes in Rust from the forest exported once from CHORD.rds (it also reads the .rds
  itself); equal to R on 20,000 random samples. 7-15 times faster, a tenth of the memory, no R. Findings
  CH1-CH10: R breaks BRCA1/BRCA2 ties at random (CH1), turns sample ids such as `0123` into `123` (CH3), and
  a soft-masked reference loses every indel's context (CH4).

- **NEO** (`neo-rs`): NEO_FINDER and NEO_SCORER byte-identical to Java on HG008, a private pair and 45 of 46
  small cases (9 with synthetic RNA inputs); the 46th is Java's thread-timing row order (NE10). 3-11 times
  faster, in a tenth of the memory. Findings NE1-NE14: wrong intron lengths in cigars (NE1, NE2), cohort TPM
  sums starting at -1 (NE6: every DNA-only run's expected TPMs are 1 too low, some negative), a missing
  score distribution read as the best rank (NE13), and IG V-segment genes not excluded (NE14).

- **CUPPA** (`cuppa-rs`): CuppaDataPrep (Java) and pycuppa's prediction (Python, scikit-learn) in one binary,
  from the classifier exported once from its pickle: byte-identical to Java and pycuppa on 30 cases (macOS arm64
  references); 10-50 times faster, a quarter of the memory or less. pycuppa's own numbers depend on the machine's
  BLAS: the Linux pipeline's differ from the macOS references in the last 2-3 digits (at most 2.4e-13 relative on
  HG008), and the port gives the macOS ones. Findings CU1-CU8: the alt-splice-junction model expects
  log-transformed counts that it never gets (CU8, RNA runs), a missing sex feature means male (CU2), and a
  sample's probabilities depend on the cohort it is predicted in (CU7). Its summary plot (R in Java) is drawn natively
  since 2026-11-27: the layout of every element within 0.01 px of R's, PDF text and word positions equal; 0.26 s
  against 4.5 s.
  The image's toolchain is pinned to rustc 1.97.1: with rustc 1.100.0, a release build of cuppa-rs mistyped its
  test VCF's SNVs (a compiler fault or undefined behaviour, not yet found).

- **qsee** (`qsee-rs`): every data file byte-identical to Java qsee 1.0 on 107 cases (gzip bytes included),
  which also equal the pipeline's qsee 1.1-beta.1 outputs on HG008 and a private pair; without the R report 0.2 s
  and 34 MB against 2.9 s and 1.36 GB. Its visual report (PNG and PDF, R in Java) is drawn natively
  since 2026-11-27, every layout element within 0.01 px of R's and the PDF's words at R's positions; 0.3 s against
  8.8 s. Findings QS1-QS12: contamination
  shown as `>` where PURPLE uses `>=` (QS4, every run), PURPLE's WARN_CONTAMINATION shown as PASS (QS5), a
  missing BQR file fatal despite `-allow_missing_input` (QS1).

- **ORANGE** (`orange-rs`): the JSON report and the plot copies byte-identical to Java ORANGE 5.1 on 33 cases (24
  identical JSONs; in 9 both fail): HG008 with the pipeline's command, a private pair, tumor-only, targeted, minimal
  and option cases, the Java tests' run with RNA, and a synthetic case reaching fusions, germline findings, viruses
  and ISOFOX results. On HG008 0.28 s and 43 MB against 9.0 s and 764 MB (1.0 s and 350 MB for Java's JSON part
  alone). Findings OR1-OR11: `-primary_tumor_doids` never resolves (OR1), CUPPA's `simpleDups32To200B` read from the
  20 kb-1 Mb deletion feature (OR5: HG008 5 instead of 1), TAPBP misspelt "TABBP" in the immune escape genes (OR9),
  and every tumor/normal run without CHORD, so every targeted tumor/normal run, failing in the PDF writer before
  the JSON is written (OR10). Since 2026-11-27 orange-rs also draws the PDF (hmf-plot-rs's PDF writer, on
  pdf-writer): in Java mode equal to Java's on all 29 cases where Java writes one, in text, every word's position,
  embedded images and rendered pixels; 0.5-1.0 s and 110-150 MB against 4-12 s and 0.5-0.9 GB. Findings OR13-OR16
  (layout faults of the PDF: values broken inside cells, titles left alone at a page's foot).

For the R and Python parts, the first version keeps calling the same scripts, so plots and
predictions stay identical. Replacing them is a separate decision per program: a random forest or a
fitted classifier can be reproduced exactly in Rust; plots and the PDF cannot be byte-identical and
would be checked by eye.

**Native plots** (`hmf-plot-rs`, from 2026-11-27): a shared drawing crate on BUI's software raster
(https://github.com/ppmpreetham/BUI, pinned by commit; its C kit is built with cmake), with the ggplot2 elements
the R scripts use and a Circos ring engine driven by the same .conf files. PURPLE's 9 plots come first: same file
names, pixel sizes and conditions as Java. In `-java_compat` the 7 R plots match Java's panels and text extents
within 1 px; the 2 Circos images differ in rendering only (anti-aliased shapes, dots as discs). The PNGs are
byte-identical across runs and across macOS arm64, Linux arm64 and Linux amd64 (libm's sin/cos/pow, C built with
`-ffp-contract=off`). The charting step takes about 0.7 s, against about 24 s for R and Circos. Findings PL1-PL4:
Circos's `skip_run` drops 15% of HG008's plotted SNVs (PL1), ggplot2 3.4+ ignores the bars' outline size (PL2),
and the rainfall distances cross chromosome boundaries (PL3). Since then LINX's visualiser, ORANGE's PDF, CUPPA's summary
plot and qsee's visual report (PNG and PDF) are drawn natively too (below and in each crate's REPORT.md), so no
figure of the pipeline needs R, Perl, Circos or Java.

### Wave C: RNA and newer tools

| Program | Size | Notes |
|---|---|---|
| ISOFOX | 27,800 | RNA transcript counts, novel splice junctions, retained introns and fusions. Ported (`isofox-rs`, 2026-11-27): byte-identical to ISOFOX 2.1 on a private tumor RNA slice and synthetic cases, 14-19 times faster in a quarter of the memory; findings IF1-IF13 (e.g. every multi-exon gene's intronic length too large, IF2). ISOFOX 2.1 expects bwa-mem3 + TARS BAMs; for STAR BAMs (Hartwig's README asks for 2.0) the port also reproduces ISOFOX 2.0.1 (`-isofox_version 2.0.1`, or chosen from a STAR @PG line), byte-identical to the 2.0.1 jar on the same cases; findings IF14-IF18 in 2.0.1-only code |
| Viridian | 4,100 | Hartwig's replacement for VIRUSBreakend and Virus Interpreter; not yet in oncoanalyser 3.0.0 |

SAGE append, used for RNA, is already part of `sage-rs`.

### Not planned

Cohort and research tools that oncoanalyser does not run: compar, patient-db, panel-builder,
pave-reverse, gene-utils, sv-tools, paddle, v-chord, v-cuppa, taur, hmf-id-generator.

Programs outside hmftools: see section 7.

## 5. Shared code

Done on 2026-11-26: `hmf-common-rs` holds what the three ports had each copied, and they use it. That
is the indexed FASTA reader and reference genome, Java number formatting, option parsing, logging,
regions and chromosome order, and the engine of the `-java_compat` switch (each port keeps its own
finding ids). It also has ESVEE's Java-behaviour helpers (HashMap iteration order, identity-hash
order under `-XX:hashCode=2`, TimSort, stream sums), which later ports need. BAM reading stays in each
port: the three read differently, each tuned to its access pattern.

Each port takes the crate as a path dependency rather than joining one Cargo workspace. That keeps
each port's `Cargo.lock`, build folder and binary path, which the test scripts use. The check was the
full set of regression runs of the three ports against the stored Java outputs, in Java mode. All
results were as before; the failures are the same ones, known in advance, from Java's own faults
or thread-count effects. In addition, the old and the new binaries gave identical output in the
default mode on one input each.

Before wave B, which leans heavily on `hmf-common` (Ensembl cache, gene panels, driver catalogues,
VCF and TSV readers), those readers go into this crate too.

## 6. How each port is checked

As for the first three:

- the Java unit tests ported;
- small slices, then the whole-genome pair, plus the public HG008 data where it applies;
- identical output to Java under `-java_compat`, compared with a Java baseline built from the same
  hmftools commit (one heavy Java run at a time on the shared machine);
- each Java fault found gets an id, is corrected by default and is listed in the port's REPORT;
- time and peak memory against Java on the same machine.

Machine time for Java baselines and whole-genome runs, not writing code, has been the slowest part of
each port so far.

## 7. Programs outside hmftools

With waves A to C done, the pipeline runs no Java, R, Python or Perl from hmftools, but it still runs
programs from elsewhere. The rule for them: a program that does a lot of the work stays (the aligners);
where the pipeline uses one or two functions of a large package, a Rust implementation of just those is
usually better. It is one small binary with no JVM or Python, it can be checked record by record against
the original, and we can adapt it.

| Program | What it does in our runs | Cost seen | Plan |
|---|---|---|---|
| bwa-mem3 | DNA alignment | most of a run's CPU time | stays |
| sambamba view and sort | sorts bwa-mem3's output for REDUX | inside the alignment tasks, not measured on its own | stays for now; measure it, then consider a sort in Rust |
| STAR 2.7.3a | RNA alignment | 14 min on 64 vCPUs for a 74 M read-pair tumor RNA | stays. Hartwig's newer RNA route (bwa-mem3, TARS, REDUX) would replace it once TARS's resources are published |
| samtools sort | sorts STAR's output | 7.5 min using 2.6 cores | to be folded into the duplicate marking below |
| GATK MarkDuplicates (Picard) | marks duplicate RNA reads before ISOFOX | 15.3 min using 1.8 cores of 64, 59 GB peak (its heap is set to 80% of the task's memory) | ported: `markdups-rs` (below) |
| samtools index; `samtools --version`, `sambamba --version` | TEAL's prep indexes its output BAMs; version lines in `versions.yml` | seconds | teal-rs to write the index itself; then the Rust image needs neither program |
| VIRUSBreakend (GRIDSS 2.13.2 with kraken2) | viral presence and integration | 1.9 CPU hours, 5.9 GB | Hartwig's Viridian replaces it and Virus Interpreter. Port Viridian when oncoanalyser adopts it. It calls BWA through JNI, and esvee-rs already vendors BWA |
| MultiQC | summary report | minutes | stays (its own container; optional) |
| gpgr's LINX report (R) | LINX HTML report | minutes | optional; not planned |
| fastp | not run in our runs (no UMIs) | none | if UMIs are needed, moving the UMI into the read name is small enough to write in Rust |

### Duplicate marking for RNA

Decision (2026-11-27): port Picard's MarkDuplicates, as GATK 4.6.1.0 runs it in oncoanalyser, to Rust.
Under `-java_compat` it gives the same duplicate flags and metrics, and it corrects Picard's faults by
default. Then fold the sort in. STAR's unsorted output keeps each pair's records together, so duplicates
can be found as the records stream in and the sorted BAM written once. That saves samtools sort and a
pass over the BAM.

Done (2026-11-28): `markdups-rs`, with the `gatk` wrapper in the image. Under `-java_compat` its output equals
GATK 4.6.1.0's (the BAM's uncompressed bytes, and the metrics apart from the date line) on a tumor RNA slice,
three public HG008 slices and 39 synthetic cases, at 1 and 2 threads; the one-step mode (`SortAndMarkDuplicates`)
equals samtools 1.21 sort followed by GATK. On the RNA slice it takes 5 s and 0.3 GB against GATK's 21 s and 2.6 GB
(2 threads). Seven Picard faults are corrected by default (MD1-MD7, in its README); one touches this kind of data:
cluster coordinates above 32767 wrap in Picard (MD1), which on the RNA slice changes 3 tied duplicate pairs and
optical duplicates from 19,851 to 19,885.

Alternatives looked at:

- **riker** (Fulcrum Genomics; Rust, MIT): successors of Picard's QC metric tools (alignment, insert size,
  whole-genome coverage, GC bias, hybrid capture, RNA-seq). It does not mark duplicates, and its outputs
  differ from Picard's on purpose. Nothing in our runs needs it. Its RNA metrics may be worth adding later,
  since oncoanalyser runs no RNA QC beyond ISOFOX's summary.
- **FastDup** (C++; Bioinformatics, 2025): Picard's algorithm, about 20 times faster than Picard by its
  authors' measurements. It finds the same duplicate sets, but the read it keeps from a set can differ
  from Picard's, and it keeps Picard's overflow of large optical-duplicate coordinates on purpose. It
  brings its own htslib and a C++ build, and its flags could not be checked against Picard's record by
  record. Its idea of working on chunks of the sorted BAM in parallel is worth borrowing.
- **fgumi `dedup`** (Fulcrum Genomics; Rust): UMI-aware duplicate marking, with a mode for reads without
  UMIs. It does not claim Picard's results, and it needs a different sort order.
- **REDUX** (already ported): Hartwig's newer RNA route marks duplicates with REDUX, after bwa-mem3 and
  TARS. ISOFOX 2.0.1's STAR route was built on Picard's flags, so REDUX there would change ISOFOX's counts.
  It fits once the RNA route moves to TARS.
