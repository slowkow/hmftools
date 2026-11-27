# COBALT in Rust: record of the port

A record of how `cobalt-rs` was checked against Java COBALT 3.1 (this repository, `mvn -pl cobalt -am
package`), the faults found in the Java code, and timings. Options and outputs are in
[README.md](README.md).

## 1. How it was checked

Java COBALT 3.1 and the port were run on the same inputs with the same options, and every output file
was compared (gzip files after decompression; `cobalt.version` by its version line, as it has the build
time). Java runs used 2 threads unless stated and `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, the
defined reference for CO3; Java's 1- and 2-thread outputs are identical, and no output differed between
Java mode and CO3 corrected on these inputs. The port ran in Java mode (`COBALT_JAVA_COMPAT=1`).

Inputs: a small private tumor/normal pair (`TN`: chr21, chr22 and parts of chr1-10),
GC_profile.1000bp.38.cnp, -ref_genome_version 38; `TO` = tumor only with -tumor_only_diploid_bed
DiploidRegions.38.bed.gz; `GO` = germline only; panels: two synthetic panel files (windows on chr21/22, random enrichments, NA values, a row off a window
start, a duplicated position; two labels); low coverage: about 5x of the small pair's chr21/22
(subsampled with samtools); the public GIAB HG008 SBX CRAMs (chr21).

"identical" = every output file identical after decompression; the gzip bytes were identical too in every
case checked (zlib level 6 as java.util.zip).

| case | options besides the inputs | Java time, peak RSS | Rust time, peak RSS | result |
|---|---|---|---|---|
| small_tn | TN | 260.3 s, 5.82 GB | 89.0 s, 0.64 GB | identical |
| small_to | TO | 129.6 s, 4.67 GB | 43.4 s, 0.48 GB | identical |
| small_go | GO | 140.1 s, 4.72 GB | 46.0 s, 0.48 GB | identical |
| small_tn_t1 | TN -threads 1 | 427.1 s, 4.68 GB | 192.9 s, 0.63 GB | identical (Java 1 vs 2 threads: identical) |
| small_tn_mq | TN chr21;22 -min_quality 40 | 102.2 s, 2.93 GB | 20.4 s, 0.13 GB | identical |
| small_tn_dup | TN chr21;22 -include_duplicates | 102.8 s, 2.18 GB | 24.1 s, 0.13 GB | identical |
| small_tn_gc | TN chr21;22 -gc_ratio_min 0.3 -gc_ratio_max 0.6 | 75.3 s, 2.17 GB | 20.6 s, 0.13 GB | identical (Java: same ratios as with the defaults, CO1) |
| small_tn_gamma | TN chr21;22 -pcf_gamma 20 | 72.9 s, 2.39 GB | 20.6 s, 0.13 GB | identical |
| small_tn_skip | TN -skip_pcf_calc | 83.9 s, 4.89 GB | 24.9 s, 0.36 GB | identical (no .pcf) |
| small_tn_regions | TN -specific_regions chr21:20000000-30000000;chr22:1-50818468 | 98.2 s, 2.15 GB | 20.6 s, 0.13 GB | identical (whole chr21 and chr22, CO6) |
| small_tn_t0 | TN chr21;22 -threads 0 | 0.4 s, exit 1 | 0.0 s, exit 1 | same failure, no files (CO5) |
| small_tn_x | TN -specific_chr chrX | 2.3 s, exit 1 | 0.5 s, exit 1 | same failure (no autosome: IllegalStateException), no files |
| small_tn_panel | TN chr21;22 -target_region_norm_file panel_one | 70.6 s, 2.20 GB | 19.9 s, 0.13 GB | identical |
| small_tn_panel2 | TN chr21;22, panel_two (2 labels) | 70.5 s, 2.43 GB | 19.9 s, 0.13 GB | identical |
| small_to_panel2 | tumor only chr21;22, panel_two | 36.5 s, 2.08 GB | 9.7 s, 0.13 GB | identical |
| small_go_panel2 | GO chr21;22, panel_two | 39.6 s, 2.14 GB | 10.4 s, 0.13 GB | identical |
| lowcov_tn | 5x TN (median 5.784: consolidation count 10) | 21.7 s, 4.61 GB | 5.0 s, 0.50 GB | identical |
| lowcov_to | 5x TO | 15.0 s, 4.78 GB | 3.3 s, 0.35 GB | identical |
| sbx_cram_tn | HG008 SBX CRAM pair, chr21, -ref_genome | 111.7 s, 3.09 GB | 7.3 s, 0.33 GB | identical |
| sbx_cram_dup | the same, -include_duplicates (consensus reads dropped) | 115.2 s, 2.93 GB | 13.5 s, 0.40 GB | identical |
| e2e (26 scenarios) | the Java ProcessBamTest scenarios on its synthetic BAMs, plus one with -ref_genome_version 37 (below) | 0.5-3 s each | < 0.1-1 s | 26 identical |

(Rust's small-pair runs shared the machine with the Java queue and other delegates' jobs; times vary by
up to 2x with load. chr21;22 = -specific_chr chr21;chr22.)

COBALT has no sequencing-type logic (no option, no UQ or SBX handling), so the HG008 data serve as the
CRAM check; the Ultima slices add nothing COBALT treats differently.

e2e scenarios (`src/e2e_tests.rs`; with COBALT_E2E_KEEP=<dir> the tests keep their synthetic inputs, on which Java
was run): tumor-only, germline-only, tumor/normal on one window (targeted), segmentation
files by mode, diploid BED filtering, the GC profile used below depth 1 (with overlapping GC profile
lines), panel enrichment with a 0 enrichment, low coverage whole genome and targeted, mappability at and
around 0.85, unmappable windows, reads spread over windows, extreme GC (targeted and whole genome),
depth-0 windows (both), GC normalisation, segmentation (targeted and whole genome), female, male, depth-0
regions in tumor and reference, the pseudo-gene exclusions of chr15/16 at full 38 lengths, and tumor/normal
with -ref_genome_version 37 (names without chr, the 37 exclusions). All 26 are identical between Java and
Rust. Java's ProcessBamTest is @Ignore'd; two of its expectations do not hold
for Java either (referenceGCDiploidRatio of windows without reads in targeted mode: Java writes -1, the
test expects 0); the port's tests check Java's actual value there.

Unit tests: `cargo test --release`: 96 pass (75 run without COBALT_TEST_RESOURCES, which skips the 21 e2e
tests; set it to `<hmftools>/cobalt/src/test/resources`, or run from `<hmftools>/cobalt-rs`, where
`../cobalt/src/test/resources` is found). Ported Java tests (`src/java_tests.rs`, one Rust test per Java
test or group): BamRatioTest, CobaltWindowTest, GCPailTest, GCPailsListTest, GcBucketStatisticsTest,
ResultsConsolidatorTest, LowCoverageConsolidatorTest, BamRatiosTest, DiploidRatioNormaliserTest (Java's
expectations hold in Java mode only: they encode CO2; the corrected values are checked too),
DiploidNormaliserTest, ReadDepthStatisticsNormaliserTest, UnityNormaliserTest, ResultsCollatorTest,
WindowStatusTest, WindowStatusesTest, SuppliedExcludedRegionsTest, ExcludedRegionsFileTest (the resources),
DiploidRegionLoaderTest, TargetRegionsTest, WholeGenomeTest, CobaltCalculatorTest (all five),
CobaltRatioSegmenterTest (with its arm locator, through the shared segmentation), TumorRatioSegmenterTest /
ReferenceRatioSegmenterTest (the value columns), ReadDepthAccumulatorTest, RollingMedianTest,
ProcessBamTest (the e2e scenarios). Not ported: tests of Java object plumbing (equals, hashCode, toString)
and of class wiring with Mockito (BamCalculationTest, TumorCalculationTest, ReferenceCalculationTest,
DoNothingNormaliserTest, NoOpConsolidatorTest), ReadCountMedianTest (hmf-common's
Doubles.interpolatedMedian, unused by COBALT), the metrics and norm tests (out of scope, section 5). Own
tests: CO2 and CO4 in both modes, the Java key order under -XX:hashCode=2 against a Java harness,
GC counting on packed bases (every range), window splitting, config parsing, file formats, panel labels.

Library checks with small Java harnesses compiled against the COBALT jar: the ArrayListMultimap key
order under -XX:hashCode=2 (three maps chained as COBALT copies them), and RollingMedian (removing 2 from
{1, 2, 3} leaves 3 values).

## 2. Findings

Corrected by default; `-java_compat` (or `-java_compat_except CO<n>,...`) gives Java's behaviour. All of
CO1-CO6 are in COBALT 3.0 too (which upstream's module files name), checked against its sources at tag cobalt-v3.0.

- **CO1. `-gc_ratio_min` / `-gc_ratio_max` do nothing.** CobaltConfig stores them in
  CobaltConstants.GC_RATIO_MIN/MAX, which nothing in COBALT reads: the GC bucket limits GC_BUCKET_MIN/MAX
  are computed from the defaults when the class loads, and ReadDepthStatisticsNormaliser uses the
  DEFAULT_ constants. (The panel builder NormalisationConfig does set GC_BUCKET_MIN/MAX from these options,
  which shows the intent.) Evidence: Java with 0.3/0.6 writes the same ratio and gc.median files as with
  the defaults (small_tn_gc vs small_tn_gamma). Corrected: the options give the bucket limits
  (calcGcBucket) and the GC range of the mean/median normaliser. Effect: none at the defaults (every
  pipeline run); the options work as described. With bucket limits of 0 or 100 the run fails as
  GcBucketStatistics' preconditions do in Java.
- **CO2. RollingMedian loses removals.** DiploidRatioNormaliser keeps the ratios within 5,000 valid
  windows in a two-heap rolling median. RollingMedian.remove takes a value from the max-heap when it is at
  most the median, else from the min-heap; PriorityQueue.remove does nothing when the value is not in that
  heap. When the min-heap holds the extra element its top is the median, so removing the median value
  searches the max-heap, finds nothing, and the value stays in the window for good (and the window keeps
  growing by one). Evidence: harness (remove 2 from {1, 2, 3}: 3 values left); Java's own
  DiploidRatioNormaliserTest asserts values computed with a stale 0.51 (the true windows give 0.49/0.535
  and 0.59/0.49 for its last two ratios); small pair: 75 referenceGCDiploidRatio values change, the .pcf
  files do not. Corrected: the value is removed from the heap that holds it. Effect: referenceGCDiploidRatio
  (which PURPLE reads) and the reference .pcf means. Whole genome (T/N, corrected vs Java mode): 152,925
  referenceGCDiploidRatio values change, by at most 0.0002; the MeanRatio of 125 of 2,650 reference
  segments changes by at most 0.0002; no segment boundary moves.
- **CO3. The mean read depth depends on identity hashes.** The read depths whose mean and median give the
  mean/median normalisation factor (and the means of UnityNormaliser in targeted mode) are summed in the
  key order of a Guava ArrayListMultimap keyed by the HumanChromosome enum, a HashMap of identity hashes:
  the order, so the last bits of the mean, can differ from run to run. Java mode reproduces the order of a
  JVM with `-XX:hashCode=2` (computed with hmf_common::java_id_hash, checked against a Java harness;
  the small-input Java references here were made that way; the whole-genome ones with the default JVM,
  and also matched); corrected: chromosome order. Effect: none seen (small pair,
  T/N: identical files either way); at most a last-digit change of a 4-decimal value in rare windows.
- **CO4. A window without reads keeps ratio 0 where its GC bucket masks every other window.**
  BamRatio.normalise returns early for a ratio of 0 whatever the factor, so a zero-depth window whose GC
  bucket has no normalisation factor (outside the GC limits, effectively buckets 25-67, or a bucket whose
  smoothed median is 0) keeps ratio 0, while every window of that bucket with reads is masked (-1). Zero
  ratios are segmented as log2 floor -9.97 and read by PURPLE as copy number 0, so this leaves spurious
  "no coverage" windows exactly where COBALT's own GC criteria say the data cannot be trusted (GC-extreme
  windows are the ones most likely to get no reads). Corrected: a ratio of 0 is kept only when the factor
  applies (positive, window included); otherwise it is masked like the others. Effect, small pair (whose
  tumor and normal have no reads on most of chromosomes 1-20, X, Y): 10,269 tumorGCRatio and 10,269
  referenceGCRatio values 0 -> -1, 3 changed segments in each .pcf. Whole genome at 88x/115x: 38 tumor and
  36 reference windows 0 -> -1 (of about 9,860 zero-ratio windows each), one tumor segment mean changes by
  0.0001, no boundary moves. It matters more at low depth, where more GC-extreme windows get no reads.
- **CO5. `-threads 0` fails.** Executors.newFixedThreadPool(0) throws IllegalArgumentException outside any
  handler: exit 1 with a stack trace, no output, although the option's description says "0 or 1 not
  multi-threaded". Corrected: taken as 1. (Same as AMBER's AM7.)
- **CO6. `-specific_regions` restricts to whole chromosomes.** BamReadCounter uses SpecificRegions only to
  drop chromosomes, so `-specific_regions chr21:20000000-30000000` processes all of chr21 (small_tn_regions:
  every window of chr21 and chr22 in the output). The description says regions. Kept as Java (no switch):
  limiting the reads to the regions would change every normalisation (medians, means, rolling windows), so
  the intent is unclear; the id is reserved for the record.

Remarks (Java kept; not faults, or intent unclear):
- GC buckets: GcBucketStatistics stores bucket i's smoothed median when bucket i + 1 is allowed, and
  allows (min, max]: with 24/68 the windows of buckets 25-67 are normalised and bucket 68 is masked
  (README says [0.24, 0.67]; Java's tests expect this). ReadDepthStatisticsNormaliser takes GC 0.24-0.68
  inclusive (fractions, not buckets).
- Germline-only runs never consolidate low-coverage windows: the reference takes the tumor's consolidator,
  which an empty tumor makes a no-op (median NaN). ReferenceCalculation has a todo about it.
- A run with no autosome window of valid ratio (e.g. `-specific_chr chrX`) fails in Java with an
  IllegalStateException (mean/median NaN); the port fails too, saying why.
- Targeted mode: the median ratio file lists all 24 chromosomes, with the median of positive ratios and
  the number of all windows (whole genome: valid windows only), so CobaltChromosomes' Y count check (1,000)
  always passes; only positive reference ratios get a referenceGCDiploidRatio (zero-depth windows -1);
  with several panels only the first panel's ratios get one, and only its statistics go to gc.median.
- Panel rows are used only at window starts (1, 1001, ...); other positions are ignored silently; a
  repeated position keeps the last row.
- DiploidRatioNormaliser's window is 5,000 valid windows on each side (by index, not by position), at
  least 1,000 of them (README: "10Mb window, minimum 1Mb readable").
- `-bam_validation` is accepted and checked for the three values but not applied (htslib reads leniently).
- README differences: "Median read depths for allosome reads" (the code uses autosomes); threads default
  4 (1), `validation_stringency` (bam_validation, default SILENT not STRICT), `use_old_segmenter` (no such
  option in 3.1; the R segmenter, RatioSegmentation, is dead code).
- The median ratio file has NaN for chromosomes without a valid window (e.g. the small pair's).

## 3. Timings and memory

Small pair and other cases: the table of section 1 (wall time, maximum resident set size; same thread count).

Where Rust's time goes (sample profiler on a 1-thread tumor-only run of the small tumor, reading only):
inflate (libdeflate) 54%, the counting 15% (before the 8-bytes-at-a-time G/C count, which cut instructions
by 14%), record copies and BGZF/htslib bookkeeping most of the rest. Single-thread reading runs at about
the speed of `samtools view -c` on the same BAM (18.6 s vs 19.2 s). GC profile 0.4 s, ratio calculation
0.5 s (chr21). The segmentation (shared code, O(n^2) per arm with a division in the inner loop that bit
equality needs) is 68 of the 89 s of the small T/N run: the small pair has ~2.5 million zero-depth windows
that are segmented like valid ones.

Whole genome (a private 88x/115x tumor/normal pair, 96 GB and 120 GB BAMs; 2 threads; Java with the
default JVM, so the references also test CO3's Java order against a JVM without `-XX:hashCode=2`):

| run | Java time, CPU, peak RSS | Rust time, CPU, peak RSS | result |
|---|---|---|---|
| tumor/normal, Java mode | 59.8 min, 5,349 s, 4.57 GB | 13.8 min, 1,362 s, 0.71 GB | identical, every file |
| tumor-only (`-tumor_only_diploid_bed`), Java mode | 32.4 min, 2,493 s, 5.05 GB | 5.4 min, 602 s, 0.54 GB | identical, every file |
| tumor/normal, corrected default | - | 13.1 min, 1,360 s, 0.71 GB | see CO2, CO4 |

A first Rust tumor/normal run that shared the disk with a Java AMBER run took 20.9 minutes.

Phases (`-log_debug` PerfStats), T/N: read counting 1,174 s in the run that shared the disk (sample
profile: the workers ~70% blocked in read()), 704 s in the next (about 0.31 GB/s of BAM at 2 threads, CPU-bound: inflate about half of it); GC profile 0.5 s;
ratio calculation 1.3 s; writing the ratio and median files 14.3 s (gzip level 6 of 3.1 million rows);
segmentation 32.8 s each for tumor and reference (2.59 million values each). Tumor-only: reading 298 s,
segmentation 31.5 s. The pipeline's Java COBALT used 4.4 CPU hours on a 100x pair (`../RUST_PORTS.md`
section 1); Rust used 23-25 CPU minutes here. Memory is the window arrays and ratios (about 0.7 GB at the
peak, in the segmentation), against 5-10 GB for Java.

## 4. COBALT 3.0 and 3.1 outputs, for PURPLE

The oncoanalyser fork runs Hartwig's `cobalt:3.1-beta.1` image, built from the tag `cobalt-v3.1-beta.1`,
whose COBALT sources are the ones this port follows. Upstream's module files name COBALT 3.0, so the
differences are recorded here.

PURPLE (CobaltData) reads `<tumor>.cobalt.ratio.tsv.gz` with CobaltRatioFile.readWithGender, which finds
all nine columns by name (chromosome, position, referenceReadDepth, tumorReadDepth, referenceGCRatio,
tumorGCRatio, referenceGCDiploidRatio, referenceGCContent, tumorGCContent; reference columns of -1 depth
are read as ratio 1), and `<sample>.cobalt.ratio.pcf` of tumor and reference with PCFFile.readPositions,
which reads the chromosome, start and end only (columns 0-2; MeanRatio is not used). It computes the
median ratios itself (MedianRatioFactory); the ratio.median and gc.median files are not read by PURPLE.

Between the cobalt-v3.0 tag and the fork's 3.1 the writers are unchanged: CobaltRatioFile,
CobaltMedianRatioFile, CobaltGcMedianFile (3.1 only dropped the reading of an older layout), SegmentsFile,
the file names, the column names and order, and the number formats. So PURPLE reads 3.1 outputs as it
reads 3.0's. Values differ: the segmentation penalty (3.0: one penalty per arm once there are 100,000
values, i.e. any whole genome; 3.1: one penalty for all arms, in arm order; the shared crate's REPORT
saw 3,610 vs 3,653 tumor segments on the test pair), targeted mode with several panel labels (new in 3.1),
and refactoring that does not change values (HumanChromosome keys). CO1-CO6 are in both.

## 5. Left out

- `norm` (NormalisationFileBuilder etc.: builds a panel normalisation file; oncoanalyser runs it only for
  that) and the `utils` classes (CobaltOutputsComparison, TargetedRegionAnalyser, RawCobaltRatio(File),
  DoubleDifference): out of scope by the brief.
- `metrics` (FragmentMetrics, PartitionReader, WindowStatisticsFile, ...): a separate entry point
  (FragmentMetrics.main) that CobaltApplication never calls.
- RatioSegmentation and `r/ratioSegmentation.R`: the old R segmentation, which nothing calls in 3.1.
- `-bam_validation` semantics (see the remarks): htslib has no ValidationStringency.
- Java's exact log text (timestamps, Double.toString, stack traces): messages and order kept, not bytes.

## 6. Open points

- CO6 is kept as Java (whole chromosomes for `-specific_regions`): restricting the reading to the
  regions would change every normalisation, and the option is for debugging, not used by the pipeline.
- The segmentation dominates small runs; a faster exact search would need care to keep Java's ties
  (hmf-common-rs REPORT).
