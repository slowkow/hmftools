# LINX in Rust: record of the port

This report covers how `linx-rs` was checked against Java LINX 2.4, the faults found in the Java code,
and the timings. Java LINX 2.4 is this repository's `linx/target/linx-2.4-jar-with-dependencies.jar`,
built from sources that follow upstream commit ca04e6c39d of 2026-11-20. Options and outputs are in
[README.md](README.md). Since 2026-11-27 the port follows upstream master `883ab11bba` (the IG fusion changes of
`7a4490caa6`): section 8 says what was taken, how it was checked and what it changes; sections 1 to 7 describe the
port as checked against `ca04e6c39d`.

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

### In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 7 on Google Batch ran LINX_SOMATIC and LINX_GERMLINE as the port, in Java mode. Every step
before them came from the first run's cache, through the same PURPLE output.
- **Against the port run locally:** all 20 files are identical to the local `-java_compat` run on the
  same inputs, which equals Java with `-XX:hashCode=2`.
- **Against the first run's Java LINX** (Hartwig's 2.4-beta.3 image, default JVM):
  - 13 files are identical: breakend, driver catalog, drivers, fusion, links, neoepitope, five of the
    visualiser files, and the germline driver catalog and links.
  - The rest differ only in LX1's ids and order: clusters, svs and vis_sv_data (cluster and chain ids),
    and the germline breakend, clusters, disruption and svs files (ids and row order). Sorted, with the
    id columns left out, they are identical.
- **Later steps:** ORANGE's JSON is unchanged. The LINX visualiser's plot of chr8, the LINXREPORT HTML
  and ORANGE's PDF differ, as they show the cluster ids.
- **Time and memory:**

  | Step | Port | Java |
  |---|---|---|
  | LINX_SOMATIC | 2.8 s, 164 MB | 6 s, 1 GB |
  | LINX_GERMLINE | 3.1 s, 89 MB | 7.5 s, 1 GB |

The rerun cost $0.25, mostly for the LINX visualiser and TEAL, which also ran from the new image.

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

- **LX6 (upstream `7a4490caa6`).** The 3' genes of the known enhancer pairs (IGH-MYC, TRB-TAL1, ...) get a larger
  pre-gene distance: 1 Mb since `7a4490caa6`, 250 kb before. `BreakendGenePrep.setSvGeneData` turns it into regions
  upstream of those genes and searches every gene from a breakend in such a region with that distance. So every gene up
  to 1 Mb from the breakend is annotated, not only the enhancer partner, and can form fusions. High-impact promiscuous
  fusions (FGFR1/2/3 5') skip the upstream-distance check, so they become reportable with 3' partners up to 1 Mb away,
  which no other rule allows (known pairs 100 kb, others 10 kb). Synthetic FGFR2 cases: 4 of 7 newly report a partner
  301-652 kb away; 3 replace the partner reported before (LMO2 at 5 kb by CSTF3 at 736 kb, NDUFS3 at 5 kb by PTPRJ at 536
  kb, SERPINB3 at 63 kb by VPS4B at 303 kb). Breakend files gain the unrelated genes (HG008 +150 rows, the private pair
  +261). The commit's own statement is "enhancer fusions set 3' gene upstream distance". Corrected: the specific distance
  for the genes it is given for, the default distance (`-fusion_gene_distance`) for the others. In a window Java also
  searches with the specific distance when `-fusion_gene_distance` is larger, missing genes between the two; the
  correction searches with the larger.

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

`cargo test --release` runs 127 tests. They are LINX's Java unit tests, ported with their synthetic
fixtures (at `883ab11bba`), and 5 tests of the port's own for what `7a4490caa6` changed (`src/tests/upstream.rs`: coding
exonic 5' breakends to a pre-gene 3' breakend and across SVs, the enhancer pre-gene distance in both modes (LX6), and
the enhancer exon regions of the visualiser files):
- clustering: SimpleMergeRuleTest, ComplexMergeRuleTest;
- chaining: ChainingSimpleTest, ChainingComplexTest, ChainRoutinesTest, LinkingTest;
- analyser: SyntheticTest, InversionPairTest, TranslocationPairTest, FoldbackTest, FilteringTest,
  CopyNumberTest, AnnotationTest, ReplicationRepairTest, LineTest, DoubleMinuteTest;
- fusion: FusionTest, FusionRulesTest, FusionTranscriptTest, SpecialFusionsTest, DisruptionTest;
- gene: GeneCollectionTest;
- drivers: AmpDriversTest, HomDisruptionsTest;
- germline: GermlineDisruptionTest;
- misc: MiscTest, JcnCalcsTest;
- visualiser: CircosConfigWriterTest, CircosDataWriterTest, ColorPickerTest, ProteinDomainColorsTest,
  ScaleContigTest, ScaleIntronsTest, DisruptedExonsTest, FusedExonsTest, SegmentsTest (with 9 tests of the
  port's own: Java's HashMap orders, COBALT centring, downsampling, R's collation).

Tests of LX2 and LX5 cover both modes.

Left out:
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
  - LINXREPORT (below).
- **No unit tests yet for LX3 and LX4.** No comparison sample reaches LX2 to LX5 (LX5, for example,
  needs an AMP driver on an acrocentric chromosome).
- **The visualiser** is ported: section 7.
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

## 7. The visualiser

`linx-rs com.hartwig.hmftools.linx.visualiser.SvVisualiser ...` runs LINX's visualiser: the plan of plots, the
Circos configuration and data files Java writes, and the images drawn natively with hmf-plot-rs (the Circos image,
`fusionPlot.R`'s fusion panels and `chromosomeRangePlot.R`'s Gviz strip, composed as the R scripts compose them).

### How it was checked

Java's visualiser (this repository's LINX 2.4 jar) ran with Circos 0.69-8 and R 4.5 from the pipeline's image
(`quay.io/biocontainers/hmftools-linx:2.3.1--hdfd78af_0`, docker, 2 CPUs, 2 threads); the port ran with
`-java_compat` at 2 threads. Data files were compared byte for byte, images by pixel statistics and side by side.
Cases: GIAB HG008 (all plots; `-plot_reportable`), HG008 with all fusions written to the vis files (13 fusion plots,
`-clusterId` and `-plot_reportable`), HG008 with empty vis files, six option sets (`-chromosome` with fragile sites
and LINE elements, `-gene`, no AMBER/COBALT, `-exact_position -show_sv_id`, `-step`, `-plot_cluster_genes`), and the
private pair (all plots; `-plot_reportable`).

- **Data files:** all 5,219 identical.
- **Images:** all 222 with Java's names and pixel sizes; per case 0.7-2.1% of pixels differ by more than 32 levels
  (anti-aliased shapes where GD draws aliased ones; text glyphs rendered by FreeType through BUI rather than by
  GD; labels within 0-2 pixels of Java's). The Circos part uses Circos's geometry, the strip Gviz's, the fusion
  panels ggplot2's (with the layout metrics of cowplot's pdf(NULL) device).
- **Determinism:** the same bytes on a second run, with 1 or 2 threads, and on macOS arm64 and Linux arm64.

| Case | Port (2 threads) | Java + Circos + R (emulated amd64, 2 CPUs) |
|---|---|---|
| HG008, all (49 plots) | 11.8 s, 0.49 GB | 582 s, 5.4 GB |
| HG008, reportable (8 plots) | 3.3 s, 0.37 GB | 101 s, 5.2 GB |
| Private pair, all (94 plots) | 23.3 s, 0.62 GB | 1210 s, 5.7 GB |

In the pipeline, Java's LINX_VISUALISER on HG008 (both calls, 8 CPUs, native) took 5-6.5 min at 10.5-12.2 GB.
Rasterising (6 s of 16.6 s at 1 thread) is serialised by BUI's global text engine; PNG encoding (6-7 s) runs in
parallel.

### Findings (corrected by default, kept with `-java_compat`)

- **LV1:** `VisProteinDomain.name()` tests the 5' UTR twice and never the 3' UTR, so 3' UTRs are drawn as a protein
  domain (own colour, legend entry "3-Prime UTR") rather than as UTR. Corrected: UTR/Non-coding.
- **LV2:** chromosomeRangePlot.R stacks the image over the strip assuming it square (`rel_heights =
  c(circosWidth, chromosomeHeight)`); with fusion panels the image is taller, so the Circos image and the fusion
  panels are shrunk to 91% and centred. Corrected: the image's own height.
- **LV3:** fusionPlot.R's protein domain legend is drawn on a black row (R 4.5.3, ggplot2 4.0.3, cowplot 1.2.0): the
  domain names are hidden. Corrected: the legend on white, its row as tall as the legend.
- **LV4:** links get the frame after their segments', and the template greys out everything past the current frame:
  in the final image the links (and connectors) of the last segment, and all links of plots without segment tracks
  (chromosome plots of simple SVs, LINE clusters), are drawn faint grey. 23 of HG008's 49 plots. Corrected: the
  last frame shows everything in colour.

Not ported: cohort vis files (`-load_cohort_files`) and `-specific_regions`.

## In the pipeline: the visualiser (GIAB HG008, 2026-11-27)

HG008 rerun 16 on Google Batch ran LINX_VISUALISER as the port (both of the step's calls, all plots and
`-plot_reportable`), on purple-rs's outputs from the cache. It wrote the 49 and 8 plots Java wrote, under the same
names, drawn natively; ORANGE then embedded the reportable ones in its PDF, whose text and layout equal Java's.
24 s and 0.9 GB, against Java's 5 min and 12 GB with Circos and R.

## 8. Following upstream to `883ab11bba` (2026-11-27)

Upstream changed LINX once between `ca04e6c39d` and master `883ab11bba`: `7a4490caa6`, "Linx IG fusions (AUS-509)". Nothing
else under `linx/`, and none of the hmf-common code LINX uses, changed. What was decided for each part of it (the
ledger is [`../UPSTREAM.md`](../UPSTREAM.md)):

| Change | Decision |
|---|---|
| Enhancer 3' gene upstream distance 250 kb → 1 Mb | accepted with a correction (LX6: the distance applies to the enhancer partner genes only) |
| `checkEnhancerFusion` moved to SpecialFusions; FusionIdGenerator | accepted (no behaviour change) |
| AllowExonSkipping, RequirePhaseMatch, RequireUpstreamBiotypes removed, with `-skip_unphased_fusions` | accepted (the option is now refused, as by Java) |
| A coding exonic 5' breakend fuses with any valid 3' transcript (it needed the same SV and an exonic or intronic 3' breakend) | accepted |
| The visualiser data draw an enhancer fusion's upstream side over the known enhancer region it was made from, for every enhancer (only IG before, over ImmuneRegions' loci) | accepted |
| Test changes | ported |

**How it was checked.**
- The regression cases of section 1 were rerun with Java at `883ab11bba`, with the same commands and `-XX:hashCode=2`:
  the four whole cases and the 14 small ones, plus `som_fusion_opts2`, which replaces the case whose option was removed.
- Synthetic cases (public Ensembl coordinates) exercise what changed:
  - IGH-MYC at 400 kb, IGH-CCND1 at 150 kb, IGL-BCL6 at 700 kb, IGK-MYC at 1.2 Mb, IGH-BCL2 past BCL2's 3' end;
  - TRA-LMO2 and BCL11B-TLX3;
  - IGH promiscuous partners, an IGH breakend just outside the locus, IGK-IGH;
  - BCR exon 14 to ABL1 on one SV, through a templated insertion, and to ABL1's promoter region;
  - FGFR2 breakends into the 1 Mb windows.
- With `-java_compat` all 36 cases are byte-identical to Java at `883ab11bba`.
- The port is 2.2-2.9 times faster than Java at one thread, in 15-30% of the memory.

**What it changes** (default mode, against the port before):
- The real cases keep every fusion, driver and visualiser file. Their breakend files gain the enhancer partners' pre-gene
  rows (HG008: CCND1, ETV6) and lose the unrelated genes the old 250 kb window had added (HG008 −17; the private pair +10
  −13).
- Without LX6 they would gain 150 and 261 rows of unrelated genes.
- In the synthetic cases:
  - IGH-MYC and IGL-BCL6 beyond 250 kb are now reported;
  - BCR-ABL1 through a templated insertion, or to ABL1's promoter region, is now reported on the canonical BCR transcript
    in frame (exon 13 to ABL1 exon 2); before, it was a non-coding BCR transcript, out of frame;
  - nothing within the old limits changes;
  - none of the distant FGFR2 fusions of LX6 appear.

**The visualiser.** It did not change upstream. Java at `883ab11bba` draws byte-identical images and data to the old
references on the 13 cases of section 7, and the port matches them as before (5,219 data files identical; 222 images).
On three synthetic enhancer cases (`@IGH`, `@IGK`, `@IGL`, `@TRA`, `@BCL11B`) the port's data files equal Java's:
528 files, 22 images with 0.5-1.5% of pixels differing by more than 32 levels.

The per-case record is in the port's working notes (`work/up_linx/port/NOTES.md` in the workspace).
