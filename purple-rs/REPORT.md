# PURPLE in Rust: record of the port

This report covers how `purple-rs` was checked against Java PURPLE 4.5, the faults found in the Java
code, and the timings. Java PURPLE 4.5 is this repository's `purple/target/purple-4.5-jar-with-dependencies.jar`,
built from sources that follow upstream commit ca04e6c39d of 2026-11-20. Options and outputs are in
[README.md](README.md).

## 1. What is ported

Every file of a tumor/normal run:
- **Fit and copy number:** purity.tsv, qc, purity.range.tsv, segment.tsv, pre_resegmented_regions.tsv,
  cnv.somatic.tsv, cnv.gene.tsv and chromosome_arm.tsv.
- **Somatic variants:** somatic.clonality.tsv and somatic.hist.tsv.
- **VCFs:** the somatic, germline, SV and germline SV VCFs.
- **Drivers:** germline_amp_del.tsv and both driver catalogs.
- **The rest:** purple.version.
- **Circos (`-circos`):** the 12 data files, the two .conf files and gaps.txt; the port runs Circos as
  Java does.

Targeted mode also writes target_region_cn.tsv. The R plots are not ported (section 6).

Not ported, as they are unused: SomaticHistogramPeaks, ResegmentationUtils, ObservedRegionFile,
GCAccumulator, ChimerismDetection and the tools package. `-pipeline_format` takes only its default,
OA_V3_1; any other value gives an error.

## 2. How it was checked

The port ran with `-java_compat` against Java's output on the same inputs, every file compared byte for
byte.
- **VCFs:** compared decompressed, and also as written: the BGZF bytes and the tabix index. Both are
  identical to Java's in every case with a defined reference below.
- **purple.version:** compared without its build date.
- **Plots:** R and Circos were stand-ins in both runs, so every data and Circos input file was written
  and nothing was drawn.

### Defined references

Two Java behaviours vary from run to run, so each needs a defined reference:
- **Kataegis ids (PU1)** depend on thread timing. The reference is Java at `-threads 1`.
- **baf.circos (PU2)** follows identity-hash order. The reference is Java with
  `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, which slows Java 3-5 times, so the timings use
  default-JVM runs.

At `-threads 2` against a default JVM, only those two outputs differ.

### Inputs

- **GIAB HG008** (public): the AMBER, COBALT, ESVEE and PAVE outputs of an oncoanalyser run of
  2026-11-26, with the pipeline's options.
- **A private whole-genome pair:** the inputs of an all-Java oncoanalyser run, copied with the sample ids
  replaced by TUMOR and NORMAL in file names and contents.
- **Germline-only** (synthetic): the HG008 tumor's AMBER, COBALT, ESVEE germline and PAVE germline
  files, copied under the normal's name.
- **Targeted:** HG008 with a BED of the driver genes' canonical coding exons, taken from the Ensembl
  cache.

### Results

Time is wall time and memory the peak resident set (`/usr/bin/time -l`), on an arm64 Mac, one run at a
time. Java ran with `-Xmx8g` and the default JVM options.

| Case | Threads | Java | Rust | Speed-up | Result |
|---|---|---|---|---|---|
| HG008 tumor/normal | 1 | 35.6 s, 2.66 GB | 11.4 s, 1.03 GB | 3.1x | all 30 files identical to the `hashCode=2` reference |
| HG008 tumor/normal | 2 | 26.4 s, 4.70 GB | 7.6 s, 1.05 GB | 3.5x | identical but KT (PU1) and baf.circos (PU2) |
| HG008 tumor-only | 1 | 28.4 s, 2.48 GB | 8.6 s, 0.71 GB | 3.3x | identical |
| HG008 germline-only | 1 | 34.2 s, 2.68 GB | 3.9 s, 0.84 GB | 8.8x | identical |
| HG008 targeted | 1 | 50.6 s, 2.71 GB | 14.3 s, 1.03 GB | 3.5x | identical, target_region_cn.tsv included |
| Private pair, tumor/normal | 1 | 50.6 s, 2.61 GB | 15.4 s, 0.84 GB | 3.3x | identical |
| Private pair, tumor/normal | 2 | 33.1 s, 2.64 GB | 9.7 s, 0.86 GB | 3.4x | identical but KT and baf.circos |

The private pair's PURPLE output from its pipeline run (Hartwig's 4.5-beta.5 image) also equals the
port's, apart from the version string, KT and baf.circos.

**Default mode against Java mode** (HG008 and the private pair, tumor/normal) changes only:
- KT, numbered in genome order (PU1);
- baf.circos, in genome order (PU2);
- germline amplification drivers, typed GERMLINE_AMP (PU9): SDHB in HG008; SDHB, AFDN and CRLF2 in the
  private pair.

**In the pipeline:** on HG008 rerun 3 every PURPLE input was byte-identical to the first run's (the
AMBER and COBALT version files aside). Even so, Java PURPLE's somatic VCF changed on 104 of 17,254
records, only in KT: a one-to-one renaming of 24 cluster ids. That is PU1 inside the pipeline.

## 3. Findings

Each is corrected by default and kept in Java mode (`-java_compat`).

- **PU1 (kataegis ids depend on thread timing).**
  - **What Java does:** with `-threads` T > 1, SomaticStream makes T SomaticVariantEnrichment tasks,
    with the chromosomes dealt out round-robin, and each task has its own KataegisEnrichment. Every
    KataegisQueue takes an id from one shared AtomicInteger when it is constructed (REV then FWD, so 2T
    ids are used before any cluster). Each cluster takes the next id when it is flushed, while the T
    tasks run in parallel.
  - **Effect:** a run's ids depend on which thread reaches its next cluster first. The clusters and
    their members are the same.
  - **Java mode:** after the 2T constructor ids, the tasks run one after the other in task order. That
    is exactly Java at `-threads 1`.
  - **Corrected:** clusters are numbered 1, 2, ... in genome order of their first variant, the same at
    any thread count.
- **PU2 (baf.circos in identity-hash order).**
  - **What Java does:** it collects the BAFs with `Lists.newArrayList(ChromosomeBafs.values())` of an
    ArrayListMultimap keyed by HumanChromosome. So the chromosome order changes from run to run.
  - **Effect:** the list is then downsampled to 25,000 points by taking every nth point, so the points
    written differ too. This was seen between the pipeline run and the reference run.
  - **Java mode:** the order under `-XX:hashCode=2`, where every identity hash is 1. All keys then share
    a bucket, and the table grows 16 -> 32 -> 64 instead of treeifying. The bin becomes a red-black tree
    at the 11th key, with every tie broken to the left, and the root moves to the front of the bin list.
    This is hmf-common-rs `java_id_hash::IdHashMap`, checked against 24 Java vectors of an
    ArrayListMultimap.
  - **Corrected:** genome order.
- **PU4 (sample names checked only with a germline run).** SampleDataFiles.hasValidSampleNames checks
  the somatic VCFs' sample names only when `runGermline()`; it should be `runTumor()`. So in tumor-only
  mode, a somatic VCF without the tumor sample is not caught. Corrected: checked whenever the tumor runs.
- **PU5 (SV leg sums in HashSet order), not reproduced.** SvLegCopyNumberChangeFactory collects the
  variants sharing a breakend position in a HashSet of StructuralVariant, whose hash includes the
  identity hashes of the VariantContexts inside. So the order in which the legs are summed per position
  varies from run to run, at most a last-bit difference. The port sums in variant order; no such case
  occurred in these samples.
- **PU6 (arm maximum starts at Double.MIN_VALUE).** ChromosomeCopyNumbers starts the arm maximum at
  Double.MIN_VALUE, the smallest positive double. So an arm whose segments are all negative reports
  maxCopyNumber 0.0000. Corrected: the start is -infinity.
- **PU8 (phased reportability tests the wrong variant).** SomaticStream.checkPhasedReportableVariants
  tests `hasPhasedEffect(variant)`, the reported variant (already known to be phased), where it means
  the neighbour, `nextVariant`. So any passing unreported neighbour in the same local phase set with a
  protein impact is marked REPORTED, phased effect or not. Corrected: the neighbour's effect is tested.
- **PU9 (germline amplifications typed as deletions).** GermlineAmpDelFinder gives germline
  amplification drivers the driver type GERMLINE_DELETION, though DriverType has GERMLINE_AMP and
  DRIVERS_PURPLE_GERMLINE lists it. Corrected: GERMLINE_AMP.
- **PU3 and PU7** were withdrawn. ImmutablePurpleQC builds its status sets with `Sets.immutableEnumSet`
  (checked in the jar's bytecode), so the QC file lists them in enum order.

Java behaviour reproduced that is not a fault, noted because it decides output order:
- `Collectors.groupingBy` fills its HashMap with `computeIfAbsent`, which puts a new key at the head of
  its bucket, where `put` appends. The somatic driver catalog follows that order
  (`JavaHashMap::compute_if_absent`).
- `Collectors.toSet()` is `new HashSet` plus `add`; it sets the germline driver gene order.

Also noted, with no effect: KataegisWindow's copy constructor sets start = end, and
ChromosomeCopyNumbers.armForSegment looks up the start position twice.

## 4. Tests

`cargo test --release` runs 61 tests, ported from 15 Java test classes with synthetic data
(`src/test_support.rs` builds records, variants, driver genes, gene and region copy numbers, COBALT
chromosomes and the default fitting config):
- somatic: KataegisQueueTest;
- germline: GermlineReportedTest, and the genotype status calls of GermlineGenotypeEnrichment;
- fitting: PurityAdjusterTest, AmbiguousBAFTest, RegionFitCalculatorTest, WeightedPloidyHistogramTest;
- region: GermlineStatusCalcsTest;
- copy number: ChromosomeCopyNumbersTest, ExtendDiploidTest, CombinedRegionTest;
- drivers: TsgDriversTest, OncoDriversTest, CopyNumberDriversTest, TsgMultiHitTest;
- SVs: breakend notation and SV types (the port's own tests).

In hmf-common-rs, the StrictMath, commons-math, Java hash order and parse_double tests compare with
vectors from Java.

Not ported yet (about 25 classes):
- segmentation: SegmentationTest, ClusterFactoryTest, PurpleSupportSegment(Factory)Test,
  SegmentRefinerTest, ExcludedRegionsRefinerTest, PCFPositionsSupplierTest;
- ResegmentationTest and ObservedRegionFactoryTest;
- amp/del and fitting: GermlineAmpDelFinderTest, FittedPurityFactoryTest, AneuploidyDetectorTest;
- inputs: AmberPointsProviderTest, CobaltDataTest;
- SomaticPurityEnrichmentTest and GeneCopyNumberBuilderTest;
- SV copy number: StructuralVariantImplied/LegsFactory/PloidyFactory and CopyNumberChangeTest;
- diploid extension: ExtendDiploidBAFTest, ExtendLongArmTest, ExtendNonDiploidTest;
- PurpleCopyNumberFactoryTest, TargetRegionsCopyNumbersTest and the regression test.

The code they cover is checked by the seven identical cases above.

## 5. Speed and memory

The port is 3.1-3.5 times faster than Java at the same thread count, and 8.8 times in germline-only
mode. It uses 2.5 to 4.5 times less memory: 0.7-1.1 GB against 2.5-4.7 GB.

Where the time goes (the private pair, `-threads 2`, 10 s):
- the purity and ploidy grid fit, about 5.6 s, most of it in commons-math's Erf.erfc (the regularized
  gamma continued fractions, kept bit-identical, so there is no cheap shortcut);
- loading, 2 s: COBALT's 3 million ratio lines (decompression and parsing), on a second thread beside
  AMBER and the VCFs;
- the germline VCF, germline SV annotation and the amp/del finder, 0.8 s.

The largest parts of memory are the COBALT ratios (3 million windows, each holding its chromosome name
as a String) and the Ensembl cache. Storing the chromosome once per list would save about a quarter.

## 6. Plots (not ported)

Java draws them in `plot/` after writing the data files. The port writes every data file the plots read,
and runs Circos when `-circos` is given, but does not run R.

| Plot | Made by | Reads |
|---|---|---|
| `<tumor>.copynumber.png` | r/copyNumberPlots.R (Rscript) | cnv.somatic.tsv |
| `<tumor>.map.png` | r/copyNumberPlots.R | cnv.somatic.tsv (minor allele copy number) |
| `<tumor>.purity.range.png` | r/copyNumberPlots.R | purity.range.tsv, purity.tsv |
| `<tumor>.segment.png` | r/copyNumberPlots.R | segment.tsv |
| `<tumor>.somatic.png` | r/somaticVariantPlots.R (when a somatic VCF is given) | somatic.hist.tsv |
| `<tumor>.somatic.clonality.png` | r/somaticVariantPlots.R | somatic.clonality.tsv, somatic.hist.tsv |
| `<tumor>.somatic.rainfall.png` | r/somaticVariantPlots.R | plot/<tumor>.somatic_data.tsv (written by Java for R and deleted afterwards; the port does not write it) |
| `<tumor>.circos.png`, `<tumor>.input.png` | the Circos binary (`-circos`) | circos/*.circos, *.conf, gaps.txt |

ORANGE reads the circos plot, so a pipeline run without Circos needs the plot drawn some other way.

## 7. Open points

- **Not ported:** the R plots, and `-pipeline_format` values other than the default.
- **PU5** is not reproduced (section 3).
- **`-java_compat` above one thread** gives one of Java's possible KT numberings: the enrichers run one
  after the other. Java itself varies from run to run (PU1).
- **Build warnings:** 20 dead-code warnings remain, from accessors and fields ported with their classes
  but not read.

## 8. Shared code

PURPLE added to `hmf-common-rs`:
- **`variant_context`:** VCF records and headers as htsjdk 2.23 holds them (VariantContext):
  - typed INFO values: read values are strings or lists, and set values keep their Java type and are
    formatted as VCFEncoder.formatVCFField does;
  - lazy genotypes, written back verbatim unless decoded, as htsjdk's LazyGenotypesContext does;
  - decoded genotypes re-encoded with GT first, keys sorted and trailing missing values dropped;
  - header lines sorted, with the contigs kept in place;
  - the BGZF stream and tabix index of `vcf::VcfFileWriter` (PAVE's), which makes the compressed bytes
    and the index Java's. The delegate's copy built the index with htslib; the merge switched it.
- **`config_builder`:** hmf-common's ConfigBuilder: items registered by the program, `hasValue` meaning
  given on the command line, path and number checks, -version and -help, the log level.
- **`java_strictmath`:** StrictMath exp, log, log10 and pow, ported from the JDK's FdLibm.java. On the
  arm64 JVM, Math.exp, log, log10 and pow equal StrictMath on 2 million random arguments each. The tests
  compare 100,000 Java results bit for bit.
- **`commons_math`:** commons-math3 3.6.1's:
  - Beta and Gamma;
  - BinomialDistribution: cumulative, probability and inverseCumulativeProbability;
  - Poisson;
  - NormalDistribution.cumulativeProbability (Erf.erfc);
  - StandardDeviation.

  It matches 12,900 Java vectors bit for bit. hmf-common-rs `poisson` (AMBER, LINX) now re-exports it.
- **`kernel_estimator`:** hmf-common's utils.KernelEstimator.
- **Additions to existing modules:**
  - `fmt::decimal_format_min`: DecimalFormat with minimum fraction digits;
  - `JavaHashMap::compute_if_absent`: a new key goes to the head of its bucket;
  - a `java_id_hash` test from Java vectors.
- **Speed-ups that keep results bit-identical:**
  - `java_text::parse_double`: an exact fast path for plain decimals of up to 15 digits (Clinger's),
    checked against the full parse on 200,000 random strings;
  - `HumanChromosome::from_name`: a byte-level parse;
  - logGamma(0.5): computed once.

Each of these is used by other ports, whose tests and regressions were rerun after the merge.
