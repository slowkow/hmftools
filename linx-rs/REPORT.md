# LINX in Rust: record of the port

This report covers how `linx-rs` was checked against Java LINX 2.4, the faults found in the Java code,
and the timings. Java LINX 2.4 is this repository's `linx/target/linx-2.4-jar-with-dependencies.jar`,
built from sources that follow upstream commit ca04e6c39d of 2026-11-20. Options and outputs are in
[README.md](README.md).

## 1. How it was checked

Java and the port ran on the same inputs with the same options, and every output file was compared byte
for byte (`linx.version` aside, as it carries the build date).
- **Java:** the local JDK 23 with `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2`. With
  constant identity hashes, the HashMaps keyed by objects iterate in a defined order (LX1), which the
  port reproduces.
- **The port:** `-java_compat`.
- **Both:** one thread, timed with `/usr/bin/time -l`.

Inputs:
- **GIAB HG008** (public), from an oncoanalyser run of 2026-11-26: PURPLE's somatic and germline SV
  VCFs and PURPLE's output folder, with the pipeline's exact LINX_SOMATIC and LINX_GERMLINE options.
- **A private whole-genome tumor/normal pair** from an all-Java oncoanalyser run, as renamed copies
  (TUMOR / NORMAL).
  - The germline VCF's sample columns have to be renamed too. If they are not, `-sample TUMOR` finds no
    genotype, and both tools write empty files that compare equal.
- **Small cases** on HG008, covering:
  - VCF subsets of a few chromosomes, somatic and germline;
  - no known fusion file, no driver panel, and no Ensembl data (both exit 1);
  - the fusion options (`-skip_unphased_fusions -fusion_gene_distance 50000 -write_all_vis_fusions`);
  - clustering options (`-proximity_distance 2000 -chaining_sv_limit 20 -no_vis_files`,
    `-proximity_distance 20000`, germline `-proximity_distance 1000`);
  - `-sample_data_dir`, with the VCF path derived from it;
  - VCFs without records, somatic and germline;
  - germline without a driver panel.

| Case | Files | Java time, peak RSS | Rust time, peak RSS | Result |
|---|---|---|---|---|
| HG008 somatic, the pipeline's options | 14 | 2.57 s, 0.77 GB | 0.80 s, 0.20 GB | identical |
| HG008 germline, the pipeline's options | 6 | 3.20 s, 0.91 GB | 1.22 s, 0.31 GB | identical |
| Private pair, somatic | 14 | 2.59 s, 0.79 GB | 0.79 s, 0.22 GB | identical |
| Private pair, germline | 6 | 3.53 s, 1.10 GB | 1.00 s, 0.27 GB | identical |
| 14 small cases | 1-14 each | 0.4-3.3 s, 0.6-1.5 GB | 0.6-1.1 s, 0.07-0.31 GB | all identical; one exits 1 in both |

The somatic files are:
- svs, clusters, links, breakend, fusion, drivers and the driver catalog;
- the six visualiser data files;
- neoepitope.

The germline files are svs, clusters, links, breakend, disruption and the driver catalog.

**The default mode** differs from Java mode only in the LX1 ids and order. Sorted, with the ids left out,
its files are identical to Java's for HG008 and the private pair, somatic and germline. LX2 to LX5 change
none of these samples' outputs; this was checked with `-java_compat_except LXn` for each.

**Default JVMs:** two Java runs with the default identity hashes agree with each other but not with the
`-XX:hashCode=2` run. The pipeline's own outputs (Hartwig's 2.4-beta.3 image, on its machine) differ from
both in the same places:
- the artifact cluster ids in clusters, svs and vis_sv_data;
- the row order and breakend ids of the germline breakend, disruption and svs files.

Sorted, with the ids left out, all of them agree (LX1).

## 2. Findings

Each is corrected by default and kept in Java mode (`-java_compat`, or `-java_compat_except` for all but
some).

- **LX1 (output order).** Java iterates HashMaps and HashSets keyed by objects without a hashCode
  (SvVarData, SvCluster) in identity-hash order. That order changes with JVM settings and between
  machines. It reaches the output in:
  - `SvFiltering.mExcludedSVs` (the artifact cluster ids of filtered SVs);
  - GermlineDisruptions' SV-to-disruptions map (germline breakend ids, disruption row order);
  - AmplificationDrivers' cluster map (drivers.tsv GAIN order);
  - the chaining maps where their order reaches the result.

  Java mode reproduces the `-XX:hashCode=2` order (`idmap.rs`); the default mode uses insertion order.
- **LX2 (GRCh38).** `CnDataLoader.expectSingleChromosome` compares the chromosome with "X" and "Y",
  which GRCh38's "chrX" and "chrY" never equal. In male samples, LOH events on X and Y are therefore kept
  on GRCh38 but dropped on GRCh37, and on GRCh38 DeletionDrivers looks for LOH rather than hom-loss
  events there. Corrected: the chr prefix is stripped before the comparison.
- **LX3.** `PseudoGeneFinder.variantMatchesPseudogeneExons` takes both gene lists from the start
  breakend (`getGenesList(true)` twice), so the end breakend's genes are never checked. Corrected: the
  end breakend's list.
- **LX4.** `ChainRuleSelector.findJcnMatchPairs` passes the SV's JCN uncertainty and the breakend's JCN
  to `jcnOverlap` in swapped order. Corrected: the arguments in the declared order.
- **LX5 (GRCh38).** `HumanChromosome.hasShortArm` compares the short names ("13", "14", ...) with the
  chromosome as given, which GRCh38's "chr13" never equals. So AmplificationDrivers can call a whole
  acrocentric chromosome gain (GAIN_CHR) on GRCh38 only. Corrected: the chr prefix is stripped.

These are kept as Java has them in both modes, as having no effect worth a finding or as choices rather
than faults:
- `ArmCluster.setFeatures` tests `(foldbackCount % 2) == 2`, which is never true; only a warning depends
  on it.
- Ported as they are:
  - `markDeletionBridge` compares the existing link with the lower SV;
  - `classifyInversionPairClusters`' overlap branch uses `lowerBe1`;
  - the ArmCluster and LINE upper-index arithmetic.
- BreakendTransData's insert-sequence adjustment never applies in gene preparation, because the insert
  sequence is set afterwards.
- `FusionFinder.transcriptPairHasFusion`: in the third alternative-phasing search only the inner loop
  breaks, so the last upstream phasing with a match sets the skipped exon counts.
- `FusionFinder.findFusions` gives fusions of known unmappable 3' pairs two ids; ids are not written.
- AmplificationDrivers' segment count is never incremented; it reaches only the cohort file.
- SpecialFusions: where the 3' gene of a promiscuous enhancer target has no canonical transcript, Java
  goes on with a null transcript and fails later. The port skips that candidate.

## 3. Tests

`cargo test --release` runs 104 tests. They are LINX's Java unit tests, ported with their synthetic
fixtures:
- clustering: SimpleMergeRuleTest, ComplexMergeRuleTest;
- chaining: ChainingSimpleTest, ChainingComplexTest, ChainRoutinesTest, LinkingTest;
- analyser: SyntheticTest, InversionPairTest, TranslocationPairTest, FoldbackTest, FilteringTest,
  CopyNumberTest, AnnotationTest, ReplicationRepairTest, LineTest, DoubleMinuteTest;
- fusion: FusionTest, FusionRulesTest, FusionTranscriptTest, SpecialFusionsTest, DisruptionTest;
- gene: GeneCollectionTest;
- drivers: AmpDriversTest, HomDisruptionsTest;
- germline: GermlineDisruptionTest;
- misc: MiscTest, JcnCalcsTest.

Tests of LX2 and LX5 cover both modes.

Left out:
- the visualiser tests (the visualiser is not ported);
- RnaFusionTest (RNA fusion matching is not ported);
- ChainingActualTest, whose fixtures (`sample_data/`) are real sample data.

Java's tests add genes and known fusion data to caches the analyser already holds. The Rust tests build
each stage's caches first and give the analyser a snapshot (`tests/fusion.rs`, FusionEnv).

The tests found two faults in the port. Both are fixed, and neither showed in any sample output:
- `ChainUtils.reverseSectionOnBreakend`'s second branch re-adds the links in reverse order; the port
  had added them in order.
- ClusterAnalyser.annotateClusters' cluster features (ClusterMetrics.populate, annotateClusterChains,
  annotateReplicationBeforeRepair) had been left out as cohort-only. But reportClusterFeatures also sets
  TotalDeleted (to the deletion bridge length when chaining gave none), and germline disruption
  reportability reads it.

## 4. Speed

The port runs 2.6 to 3.5 times faster than Java at one thread, in a quarter to a third of the memory.
LINX uses threads only to split several samples, so one sample runs on one thread in both.

`-perf_debug` logs each phase:
- **HG008 somatic, about 1 s in all:**
  - the Ensembl data, 0.76 s, with the protein domain and splice acceptor files, as fusions are called;
  - the SV VCF 0.004 s and the breakend genes 0.004 s;
  - clustering and chaining 0.001 s;
  - transcripts, cluster annotation and drivers 0.03 s;
  - fusions 0.007 s and output 0.002 s.
- **HG008 germline (17,177 SVs), 1.1 s:**
  - the Ensembl data 0.58 s (canonical transcripts);
  - the VCF 0.25 s and the breakend genes 0.13 s;
  - clustering and chaining 0.07 s;
  - the rest 0.06 s.

Most of the remaining time is reading the Ensembl files (shared code in hmf-common-rs, which a parallel
reader could shorten).

## 5. Left out, and open points

- **Not ported:**
  - several samples per run, sample id files, and the cohort and batch outputs (`-write_cohort`);
  - RNA fusion matching (`-rna_fusions_file` is refused);
  - the diagnostic logging (chain diagnostics, `-log_verbose`);
  - the visualiser and LINXREPORT (below).
- **No unit tests yet for LX3 and LX4.** No comparison sample reaches LX2 to LX5 (LX5, for example,
  needs an AMP driver on an acrocentric chromosome).
- **The visualiser** (`com.hartwig.hmftools.linx.visualiser.SvVisualiser`, in the same jar) reads:
  - a sample's vis_sv_data, vis_segments, vis_copy_number, vis_gene_exon, vis_protein_domain and
    vis_fusion files (all written by the port), or the cohort LNX_VIS_*.tsv files;
  - the Ensembl data;
  - when given, the AMBER BAF and COBALT ratio files.

  Per cluster or chromosome it:
  - writes Circos configuration and data files: `<id>.circos.NNN.conf`, with its karyotype, link,
    segment and histogram text files;
  - runs Circos for `<id>.NNN.png`;
  - runs R (`r/fusionPlot.R`, `r/chromosomeRangePlot.R`) for the fusion and chromosome range images.

  In oncoanalyser it stays on the Java image, reading the port's files, until it is drawn natively.
- **LINXREPORT** (an R package) reads the same vis files plus the svs, clusters, links, breakend,
  fusion, drivers and driver catalog files, and writes an HTML report.

## 6. Shared code

LINX added to `hmf-common-rs`:
- `sv`:
  - StructuralVariantType;
  - a VCF record reader with htsjdk's attribute rules;
  - StructuralVariantFactory and EnrichedStructuralVariantFactory;
  - StructuralVariantData (convertSvData).
- `purple_files`:
  - PurpleCopyNumber;
  - PurityContext (purity, ploidy, gender, WGD);
  - GeneCopyNumber;
  - DriverCatalog (read and write);
  - `write_lines` (Files.write).
- `poisson`: commons-math3's PoissonDistribution, moved from amber-rs (which now uses it from here).
- `fmt::decimal_format_fixed`: DecimalFormat("0.0000"), with fixed fraction digits.
- `transcript_utils`: coding bases and phases (TranscriptUtils).
- `known_fusion`: KnownFusionCache, KnownFusionData, KnownFusionType and FusionCommon's constants.
- `ensembl`:
  - protein domain and splice acceptor data;
  - findGeneRegions, findGenesByRegion, getProteinDomainPositions and getExonRankings;
  - the biotype constants;
  - in-memory caches for tests.
