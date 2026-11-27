# ORANGE in Rust: record of the port

This report covers how `orange-rs` was checked against Java ORANGE 5.1, the faults found in the Java code, and the
timings. Finding ids: OR. Options, outputs and tests are in [README.md](README.md).

`orange-rs` writes ORANGE's JSON report (`<tumor>.orange.json`) and copies the plots into `plot/`, as the jar does.
Under `-java_compat` the JSON was byte-identical to the jar's in every case run (33 comparisons: 24 identical JSONs,
and 9 runs where both exit 1). The PDF is not written; section "PDF (not ported)" lists what Java draws and from what.

## How it was checked

Reference: `orange/target/orange-5.1-jar-with-dependencies.jar`, built from this repository's sources (upstream
ca04e6c39d of 2026-11-20).

Each case ran with the jar (`java -Xmx8g -jar`, default JVM) and with `orange-rs -java_compat`, in the case folder
(the JSON's variant and breakend plot paths are the configured directory strings), one run at a time. The JSON was
compared byte for byte, and the copied `plot/` folders file by file. Where the jar fails in its PDF writer (OR10,
OR11) the JSON-only harness gave the reference (`<case>:json`, compared with `-java_compat_except OR10,OR11`).
`-sampling_date` is given in every case but `hg008_location` (an invalid date, so today in both; see "Run-dependent
values").

Time (wall, s) and peak RSS (MB), `/usr/bin/time -l`, arm64 Mac. Java's time is mostly the PDF: the JSON-only harness
takes about 1 s.

| Case | What | Java | Rust | Result |
|---|---|---|---|---|
| hg008 | HG008 tumor/normal, the pipeline's command (`commands/ORANGE.command.sh`) | 9.03 s, 764 MB | 0.28 s, 43 MB | identical, 53,911 bytes, 21 plots |
| hg008 :json | the same, JSON-only harness | 1.00 s, 350 MB | 0.27 s, 38 MB | identical |
| hg008_tumor_only | no reference (no germline, CHORD, CUPPA, PEACH) | 7.53 s, 697 MB | 0.22 s, 32 MB | identical |
| hg008_panel | `-experiment_type PANEL`, tumor/normal | 1.75 s, 374 MB | 0.27 s, 38 MB | both exit 1 (OR10) |
| hg008_panel :json | the same, harness | 1.01 s, 350 MB | 0.27 s, 38 MB | identical |
| hg008_panel_tumor_only | PANEL, tumor-only | 7.67 s, 696 MB | 0.23 s, 33 MB | identical |
| hg008_minimal | only the required options: no LILAC, CHORD, CUPPA, PEACH, SIGS, virus, QSee, plots of SAGE or LINX, version file | 4.17 s, 533 MB | 0.22 s, 33 MB | identical |
| hg008_minimal_reference | minimal with a reference (no CHORD) | 1.82 s, 372 MB | 0.28 s, 46 MB | both exit 1 (OR10) |
| hg008_minimal_reference :json | the same, harness | 1.08 s, 347 MB | 0.28 s, 43 MB | identical |
| hg008_options | `-sampling_date 240230` (SMART resolution: 2024-02-29), `-primary_tumor_location`, `-panel_name`, `-display_sample`, `-output_id`, `-sequencing_type SBX`, `-rna_sample` (not in the VCF, no ISOFOX), `-log_debug` | 9.44 s, 750 MB | 0.29 s, 38 MB | identical |
| hg008_location | invalid `-sampling_date` (warning, today), primary tumor location | 9.40 s, 770 MB | 0.29 s, 37 MB | identical |
| hg008_doids | `-primary_tumor_doids` with the DOID JSON (OR1) and a location | 9.10 s, 755 MB | 0.28 s, 37 MB | identical |
| hg008_root | `-pipeline_sample_root_dir` (OA_V3_1 layout, defaults for every directory) | 6.63 s, 682 MB | 0.28 s, 42 MB | identical |
| hg008_root_format | root directory with `-pipeline_format OA_V2_0` and explicit plot directories | 8.31 s, 747 MB | 0.28 s, 37 MB | identical |
| hg008_sample_data_dir | `-sample_data_dir` (every file in one folder) | 9.14 s, 755 MB | 0.28 s, 40 MB | identical |
| hg008_format_file | root directory with `-pipeline_format_file` (a file naming 10 of the tool directories) | 6.22 s, 673 MB | 0.27 s, 37 MB | identical |
| hg008_format_file, first version | the same file without the CUPPA and PEACH entries | 0.98 s, 349 MB | 0.27 s, 37 MB | both exit 1 (the CUPPA file looked for in the root) |
| hg008_bad_format, hg008_bad_type, hg008_no_tumor | invalid `-pipeline_format`, `-experiment_type`, no `-tumor` | 0.3 s, 92 MB | 0.00 s, 2 MB | both exit 1 |
| hg008_no_plots | a PURPLE plot directory without the plots | 0.84 s, 306 MB | 0.20 s, 32 MB | both exit 1 |
| hg008_wrong_reference | a reference id not in the VCFs | 0.43 s, 149 MB | 0.00 s, 3 MB | both exit 1 |
| t0 | a private whole-genome tumor/normal pair (renamed copies), the pipeline's command | 7.29 s, 653 MB | 0.17 s, 36 MB | identical, plots too |
| t0 :json | the same, harness | 0.82 s, 272 MB | 0.17 s, 45 MB | identical |
| t0_tumor_only | the private pair, tumor-only | 6.04 s, 514 MB | 0.12 s, 32 MB | identical |
| testrun_wgts | the Java tests' synthetic run (GRCh37, plain VCFs, ISOFOX RNA, old LINX fusion format, DOIDs with the example DOID file) | 6.16 s, 914 MB | 0.53 s, 35 MB | identical |
| testrun_wgts :json | the same, harness | 1.12 s, 606 MB | 0.55 s, 34 MB | identical |
| testrun_panel | its targeted tumor-only run | 4.96 s, 896 MB | 0.48 s, 29 MB | identical |
| testrun_rna_tumor_id | ISOFOX without `-rna_sample` (the tumor id's files, by fallback) | 1.83 s, 589 MB | 0.53 s, 34 MB | both exit 1 (OR10) |
| testrun_rna_tumor_id :json | the same, harness | 1.11 s, 509 MB | 0.54 s, 34 MB | identical |
| synthetic | the synthetic case with RNA | 2.60 s, 430 MB | 0.26 s, 42 MB | both exit 1 (OR11) |
| synthetic :json | the same, harness | 0.97 s, 353 MB | 0.26 s, 37 MB | identical, 67,767 bytes |
| synthetic_tumor_only :json | tumor-only | 0.87 s, 288 MB | 0.21 s, 32 MB | identical |
| synthetic_no_rna | without ISOFOX | 10.96 s, 881 MB | 0.26 s, 42 MB | identical |

The pipeline's own HG008 JSON (Hartwig's ORANGE image, a different build) equals the jar's but for the key
`viridian` (the jar: `virusInterpreter`) and the order of the LINX driver plots (OR3).

**Default mode** (every case without `-java_compat`) against Java differs only where a correction applies: the plot paths (OR2) in
every case; the LINX plot order (OR3) where the directory order is not sorted (HG008); `simpleDups32To200B` (OR5:
HG008 5 to 1, the private pair and the test run too); the DOIDs (OR1); `rnaId` (OR6); in the synthetic case the
other transcripts' exon, codon and coding effect (OR7), the fusions' junction copy numbers (OR8) and the SAGE plot of
a variant whose name is a prefix of another's (OR4). OR9 changes none of these samples (tested).

**Speed.** Java's run is mostly the PDF (PDFBox drawing the large PNGs): without it, 0.8-1.1 s. The port takes
0.12-0.55 s and 29-46 MB, 3-5 times faster than Java's JSON part and 20-40 times faster than the whole jar, in a
tenth of the memory or less. Of the port's time, 85-93% is reading PURPLE's somatic VCF (HG008: 17,254 records,
0.22 s of 0.28 s; decompression is 0.02 s, the rest the parse of each PASS record's INFO and genotypes); LINX's files
take 0.03 s, everything else (interpretation, plot copies, JSON) about 0.01 s. Filtered records are skipped before
parsing in the PASS-only somatic read (as sigs-rs does); a malformed filtered record would fail in Java and not here.

## Run-dependent values in the JSON, and how Java is made defined

- `samplingDate`: `LocalDate.now()` (system time zone) unless `-sampling_date yyMMdd` is given; the pipeline does not
  pass it. Defined reference: `-sampling_date` (used in the cases), or both tools on the same day (the port takes the
  local date the same way). An unparsable `-sampling_date` falls back to today, with a warning, in both.
- Plot paths: relative to `-output_dir`, so independent of it; the variant and breakend `plotFilename`s are the
  configured `-sage_plot_dir` / `-linx_plot_dir` strings plus the file name, so they depend on how the directories
  are named on the command line. The case folders reproduce the pipeline's names.
- Order of `linxDriverPlots`, and the choice among several plots matching a prefix: the file system's directory order
  (OR3). Java and the port get the same `readdir` order on the same file system.
- Nothing else in the JSON depends on the run: no HashMap of objects reaches it (the phase set ids are a HashSet of
  Integer, in Java's order), no thread timing.

## Findings (OR)

Corrected by default, Java's behaviour with `-java_compat` (unless said otherwise).

- OR1: configured primary tumor DOIDs are never resolved. `formConfiguredPrimaryTumorDoid` looks each DOID up in a
  set that is never filled (the `-doid_json` file is not read at all), so every `-primary_tumor_doids` value logs
  "Could not resolve doid" and `configuredPrimaryTumor` stays empty; `-primary_tumor_location` is then ignored too
  (else-if). Corrected: the DOIDs are looked up in the DOID JSON (the nodes of the doid.owl graph, as
  DiseaseOntology reads them), each resolved one written once, in the order of the configured set.
- OR2: plot paths start with "/" (`"/plot/HG008-T.circos.png"`): FileBasedPlotManager joins `outputDir + "/" + "plot"`
  with an output dir that already ends in "/", and `relativePath` strips `outputDir` from `out//plot/x.png`, leaving
  `/plot/x.png`, which reads as an absolute path. ORANGE's own PDF still finds the files (PlotPathResolver falls back to
  `outputDir + "/" + path`). Corrected: `plot/x.png`.
- OR3: the LINX driver plots (and the SAGE variant plots) are listed with `File.list()`, whose order is the file
  system's: the pipeline (Linux) and this Mac give the five HG008 LINX plots in different orders. Corrected: sorted by
  file name. Java mode: directory order (`read_dir`, the same `readdir` order Java gets).
- OR4: plots are looked up by containment of a name prefix: `findReportableLinxPlot` takes the first plot whose name
  contains ".cluster-<id>", so cluster 4 gets cluster 40's (or 41's ...) plot when that is listed first, in the
  breakends', fusions' and gains/deletions' plotFilename; likewise a SAGE variant C>A gets the page of a variant C>AT
  at the same position (`findVariantPlot`). Corrected: the cluster id, or the variant, must end there.
- OR5: CUPPA's `simpleDups32To200B` is read from the feature `sv.SIMPLE_DEL_20KB_1MB` (deletions of 20 kb to 1 Mb)
  instead of `sv.SIMPLE_DUP_32B_200B`, which the vis_data file has: HG008 gets 5 where the duplication count is 1.
  Corrected: `sv.SIMPLE_DUP_32B_200B`.
- OR6: `rnaId` is never set (the builder call is missing), so it is null even with `-rna_sample`. Corrected: the
  configured RNA sample id.
- OR7: a variant's other reportable transcripts (`otherImpacts`) get the canonical transcript's affected codon, exon,
  splice region and coding effect (createOtherImpact), although the reportable effect carries its own coding effect
  and PAVE_TI its own exon and codon. Corrected: the transcript's own (PAVE_TI, and the reportable effect's coding
  effect); where PAVE_TI has no entry for the transcript, the canonical values as Java.
- OR8: a fusion's `junctionCopyNumber` is the mean of its two breakends' undisrupted copy numbers (LINX 2 dropped the
  fusion's own JCN column and ORANGE substituted these). Corrected: the junction copy number of the fusion's SV,
  0.5 x (JCN min + JCN max), as ORANGE computes it for the breakends.
- OR9: the antigen presentation genes of the immune escape analysis list "TABBP" (no such gene) for TAPBP, so TAPBP is
  never checked. Corrected: TAPBP.
- OR10: every tumor/normal run without a CHORD result fails: the PDF front page calls `hrStatus()` on the null CHORD
  record (FrontPageTables.hrDeficiencyString). That is every targeted (`-experiment_type PANEL`) tumor/normal run
  (CHORD is not loaded in targeted mode) and every run without `-chord_dir`. As the PDF is written before the JSON, no
  JSON is written either. The port writes no PDF; in Java mode it fails at the same point (exit 1).
- OR11: with RNA, a run where one reported DNA fusion has RNA support and another has none fails: the PDF's fusion
  table adds the RNA column once any fusion has support and then reads `rnaSupport()` of every fusion
  (DnaFusionTable.build). No JSON is written. Java mode fails the same way.

Not corrected (a data model change, or the intent is unclear), as Java:
- A germline aberration TERTRASOMY_9 in PURPLE's QC (hmf-common has it, the data model's PurpleGermlineAberration
  does not), and germline amp/del statuses LIKELY_DIPLOID, EXCLUDED and CENTROMETIC (not in PurpleGermlineStatus), end
  the run with "No enum constant"; so does an Isofox splice junction of type CANONICAL. The port fails the same way.
- Without a pipeline version file `pipelineVersion` is "OA_V3_1", the name of the default pipeline directory layout.
- The SAGE plot directory's default (with `-pipeline_sample_root_dir`) is under LINX's somatic directory.
- A tool missing from a `-pipeline_format_file` gets the directory "" and so the root directory itself: ORANGE then
  reads e.g. `<root>/<tumor>.cuppa.vis_data.tsv` and fails, or warns "invalid RNA sample ID" for ISOFOX.
- Every somatic deletion driver has `geneRange` "PARTIAL"; only AMP drivers are "FULL".
- The arm aberrations' `copyNumber` is the arm's median copy number, while their type and `relativeCopyNumber` use
  the mean.
- Without `-qsee_dir` the QSee plot is looked for at "null<tumor>.qsee.vis.report.png" (relative to the working
  directory); likewise the CUPPA plot without `-cuppa_dir` (with a warning).
- `localPhaseSets` lists the phase sets a variant shares with other reported variants within 1,000 bases only.
- `lowExpressionGenes` is never filled; `highExpressionGenes` are Isofox's reported genes with cohort data.
- An allele without a LILAC QC line for its gene class, a breakend outside every cytoband, two SV annotations with
  one svId, a PEACH allele count other than 1 or 2, or a CUPPA cancer type without the main classifier's probability
  end the run (NullPointerException, IllegalStateException, IllegalArgumentException), in both.
- Other failures inside Java's PDF writer than OR10 and OR11 are not reproduced (none was met).

## Tests

- 34 tests in `orange-rs` (`cargo test --release`): the Java unit tests of the ported classes, and tests of the port.
  Several of Java's immune escape and Isofox tests are `@Ignore`d and out of date: the active one is ported, and the
  current rules have tests of their own.
- `hmf-common-rs`: tests for each module added (below).

## Shared code: additions to hmf-common-rs

All additions; no existing behaviour changed.

- `linx_files` (new): LinxSvAnnotation, LinxBreakend, LinxCluster, LinxFusion, LinxDriver and LinxGermlineDisruption
  readers and file names; the breakend coordinate helpers; `Row` (a TSV line with Java's field access and errors) and
  `enum_value`, used by the other new readers.
- `small_variant` (new): SmallVariantFactory (PASS filter, HumanChromosomeFilter, NTFilter, the AD and depth checks,
  the sample checks), SmallVariant with VariantContextDecorator's fields, VariantImpact, VariantTranscriptImpact,
  AltTranscriptReportableInfo, VariantEffect, PathogenicSummaryFactory, and determineCodingEffect / CodingEffect.effect.
  purple-rs has its own copies of VariantImpact and the pathogenic summary in `purple-rs/src/variant.rs`, which could
  move to this module.
- `cytobands` (new), with the resources `resources/refgenome/cytoBands.{37,38}.tsv` (copied from hmf-common).
- `tool_results` (new): CHORD's prediction, PEACH's best haplotypes (Java and old Python formats), SIGS' allocations
  and the signature etiologies (resource `resources/sigs/signatures_etiology.tsv`), Virus Interpreter's annotated
  viruses, CUPPA's vis_data, PipelineVersionFile.
- `purple_files`: `read_purity_context_with_qc` (the PurplePurity reader qsee-rs added is shared), ChrArmCopyNumbersFile,
  GermlineAmpDel, `human_chromosome` (HumanChromosome.fromString), the germline driver catalog file name.
- `hla`: LilacQcData (`read_lilac_qc`, `lilac_qc_filename`).
- `rna`: RnaStatisticFile, GeneExpressionFile, RnaFusionFile and NovelSpliceJunctionFile readers.
- `lib.rs`: the new modules and their line in the crate documentation.

## PDF (not ported): what Java draws

`orange/report` (PDFBox with boxable; about 5,500 lines). A4 pages; each page has a header (the "ORANGE Report" title
in orange shades and the ORANGE circos logo, resource `orange_circos.png`), a side panel and a footer with page numbers
(second pass) and, with `-add_disclaimer`, a disclaimer; `-display_sample` replaces the tumor id there. Chapters, in
order (from the JSON record, `plots.*` resolved against the output directory):

1. **Front Page**: Sample summary (Samples, Primary Tumor, Date Analysed, Pipeline Version, Genome Version, Sequencing
   Type) and Technical summary (QC, Fit Method, Purity, Ploidy); Driver Summary (somatic variant, copy number and
   disruption drivers, germline variant, copy number and disruption drivers, fusion drivers, viral presence, whole
   genome duplicated, DPYD status) beside Genome Wide Biomarkers (MSI indels per Mb, TMB per Mb, TML, HR deficiency
   score from CHORD (OR10), LOH proportion, number of SVs, CUPPA cancer type); the final PURPLE circos plot
   (`purpleFinalCircosPlot`).
2. **Somatic Findings**: Small Variants (SomaticVariantTable: variant, HGVS, copies, AF, clonal, somatic likelihood,
   biallelic, hotspot, driver, RNA depth when present, the SAGE vis page names), Amplifications and Deletions
   (GainDeletionTable: location, gene, range, min/max/arm/relative CN, TPM, percentile and fold change with RNA),
   Fusions (DnaFusionTable: fusion, junctions, JCN, phasing, type, RNA support (OR11), driver), Disruptions
   (DisruptionTable: location, gene, range, type, JCN, undisrupted CN), Viruses (reportable; not tumor-only), Arm Copy
   Number Aberrations (ChrArmCopyNumberTable), Signature Allocations (not tumor-only), and the Structural Driver Plots
   (`linxDriverPlots`, two per row with their titles).
3. **Germline Findings** (not tumor-only): Small Variants (GermlineVariantTable: with ClinVar and gnomAD), Amplifications
   and Deletions, Disruptions, Chromosomal Aberrations (the germline aberrations), Pharmacogenetics (PEACH).
4. **Immunology**: HLA Class I Alleles and, for whole genomes, Class II (HLAAlleleTable: allele, QC, tumor CN,
   reference/tumor/RNA fragments, somatic mutation counts), Genetic Immune Escape (ImmuneEscapeTable).
5. **RNA Findings** (with ISOFOX): QC (RnaStatisticsTable), High Gene Expression (ExpressionTable), Fusions Detected In
   RNA, Not In DNA (RnaFusionTable), Novel Splice Junctions (NovelSpliceJunctionTable).
6. **Tissue of Origin** (not tumor-only): the CUPPA summary plot (`cuppaSummaryPlot`), or "not available".
7. **Quality Control** (when `qSeePlot` exists): the QSee report plot.
8. **Purity and Ploidy**: six PURPLE plots, landscape, 3 x 2: input circos, copy number, somatic clonality, purity
   range, minor allele map, somatic rainfall (`purple*Plot`). `purpleVariantCopyNumberPlot` (somatic.png) is copied
   but not drawn.

Data used: only the JSON record and the copied plots (`plot/`), plus the config's display id and disclaimer flag. The
QSee report PDF (`qsee/<tumor>.qsee.vis.report.pdf`) is not used; its PNG is.

## Open points

- Upstream's change of 2026-11-23 (f5f660cc9f): the virus section is `viridian`, and ORANGE also reads Viridian's
  files. The pipeline's ORANGE image already has it.
- The PDF: to be drawn natively from the JSON and the copied plots, with the shared drawing crate the plot ports use.
- In the pipeline, oncoanalyser runs ORANGE only when PURPLE's `plot/<tumor>.circos.png` exists. Until purple-rs draws
  its plots, a run with purple-rs skips ORANGE, and orange-rs runs only beside the Java PURPLE.
- `small_variant` in hmf-common-rs repeats code that purple-rs has in `src/variant.rs` (VariantImpact, the pathogenic
  summary); purple-rs could use the shared module.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 12 on Google Batch ran ORANGE as the port, in Java mode, on the inputs of rerun 11 (every other step
from the cache; PURPLE in Java, which draws the plots ORANGE needs). The JSON equals the one Hartwig's ORANGE image
(5.1-beta.1) wrote from the same inputs but for one key: the image names the virus section `viridian`, while the
5.1 jar this port follows (hmftools ca04e6c39d) and the port name it `virusInterpreter`. Upstream renamed it on
2026-11-23 (f5f660cc9f, "support loading Viridian and VirusInterpreter files"); following that is an open point.
The LINX driver plots came in the same order (in Java mode both take Linux's directory order). ORANGE took 5.3 s
and 42 MB, against Java's 18 s and 419 MB.
