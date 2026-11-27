# Rust ports of hmftools: plan for the rest of the suite

Goal: run the whole oncoanalyser chain with Rust programs, so that more samples can be processed per
machine and per day. Each port takes the same options and writes the same files as the Java program,
corrects the faults found in Java by default, and reproduces Java exactly under `-java_compat`.

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
  deletions in Java (PU9). Findings PU1-PU9. The R plots are not drawn, and oncoanalyser skips ORANGE without
  them, so the pipeline runs PURPLE as a port only with its own profile for now.
- **LINX** (`linx-rs`): LINX_SOMATIC and LINX_GERMLINE byte-identical to Java on HG008 and on a private
  whole-genome pair (14 and 6 files) and on 14 small and option cases; 2.6-3.5 times faster at one
  thread, in a quarter to a third of the memory (0.2-0.3 GB against 0.8-1.1 GB). Java's output order
  depends on identity hashes (LX1); two of the faults are chromosome-name checks that never match on
  GRCh38 (LX2, LX5). Findings LX1-LX5. The visualiser (Circos and R) and LINXREPORT are not ported.
  Inside the pipeline on HG008 (rerun 7): identical to the port run locally, and to the Java run apart
  from LX1's ids and order; ORANGE's JSON unchanged.

For the R and Python parts, the first version keeps calling the same scripts, so plots and
predictions stay identical. Replacing them is a separate decision per program: a random forest or a
fitted classifier can be reproduced exactly in Rust; plots and the PDF cannot be byte-identical and
would be checked by eye.

### Wave C: RNA and newer tools

| Program | Size | Notes |
|---|---|---|
| ISOFOX | 27,800 | RNA transcript counts and fusions; not in the run above, so its cost is not measured yet |
| Viridian | 4,100 | Hartwig's replacement for VIRUSBreakend and Virus Interpreter; not yet in oncoanalyser 3.0.0 |

SAGE append, used for RNA, is already part of `sage-rs`.

### Not planned

Cohort and research tools that oncoanalyser does not run: compar, patient-db, panel-builder,
pave-reverse, gene-utils, sv-tools, paddle, v-chord, v-cuppa, taur, hmf-id-generator.

Programs outside hmftools that stay: bwa-mem3 with sambamba (alignment and sorting), samtools,
VIRUSBreakend (GRIDSS with kraken2, until Viridian replaces it) and STAR for RNA. With waves A to C
done, the pipeline needs no Java, R, Python or Perl, but it still runs these programs.

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
