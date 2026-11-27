# ORANGE in Rust: record of the port

This report covers how `orange-rs` was checked against Java ORANGE 5.1, the faults found in the Java code, and the
timings. Finding ids: OR. Options, outputs and tests are in [README.md](README.md).

`orange-rs` writes ORANGE's PDF report (`<tumor>.orange.pdf`) and JSON report (`<tumor>.orange.json`) and copies the
plots into `plot/`, as the jar does. Under `-java_compat` the JSON was byte-identical to the jar's in every case run
(52 comparisons: 38 identical JSONs, and 14 runs where both exit 1, at the same point), and the PDF equal to the jar's
in every case where Java writes one (29 cases and a synthetic layout: the same pages, text, word positions, images and
page pixels); section "PDF" says how it is drawn and checked.

## How it was checked

Reference: the ORANGE 5.1 jar-with-dependencies built from upstream master 883ab11bba (2026-11-27), which has
f5f660cc9f ("Orange - support loading Viridian and VirusInterpreter files", 2026-11-23). The port first followed
ca04e6c39d (2026-11-20); between the two, f5f660cc9f is the only change to ORANGE's JSON or options (the other
commits touching `orange`, `orange-datamodel` or `hmf-common` only rename the modules in their pom.xml; the
hmf-common changes ORANGE uses are those of f5f660cc9f: PipelineToolDirectories(File) and
VirusInterpreterDataLoader). The 33 comparisons of the first version were all rerun against this jar, and 19 cases
were added for the Viridian input (NEW below).

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
| hg008 | HG008 tumor/normal, the pipeline's command (`commands/ORANGE.command.sh`, `-virus_dir virusinterpreter`) | 9.40 s, 752 MB | 0.28 s, 42 MB | identical, 53,903 bytes |
| hg008 :json | the same, JSON-only harness | 1.06 s, 350 MB | 0.28 s, 42 MB | identical, 53,903 bytes |
| hg008_tumor_only | no reference | 7.99 s, 696 MB | 0.23 s, 33 MB | identical, 37,596 bytes |
| hg008_panel | `-experiment_type PANEL`, tumor/normal | 1.80 s, 380 MB | 0.27 s, 37 MB | both exit 1 (OR10) |
| hg008_panel :json | the same, harness | 1.00 s, 350 MB | 0.28 s, 38 MB | identical, 39,362 bytes |
| hg008_panel_tumor_only | PANEL, tumor-only | 7.95 s, 703 MB | 0.21 s, 32 MB | identical, 37,592 bytes |
| hg008_minimal | only the required options (no virus directory: `viridian` null) | 4.01 s, 546 MB | 0.22 s, 33 MB | identical, 28,469 bytes |
| hg008_minimal_reference | minimal with a reference (no CHORD) | 1.77 s, 376 MB | 0.27 s, 43 MB | both exit 1 (OR10) |
| hg008_minimal_reference :json | the same, harness | 1.02 s, 346 MB | 0.28 s, 37 MB | identical, 29,772 bytes |
| hg008_options | `-sampling_date 240230`, location, panel, display id, `-output_id`, SBX, `-rna_sample`, `-log_debug` | 9.17 s, 757 MB | 0.28 s, 42 MB | identical, 53,919 bytes |
| hg008_location | invalid `-sampling_date` (today), location | 9.19 s, 749 MB | 0.28 s, 40 MB | identical, 53,961 bytes |
| hg008_doids | `-primary_tumor_doids` (OR1) and a location | 9.03 s, 749 MB | 0.28 s, 40 MB | identical, 53,903 bytes |
| hg008_root | `-pipeline_sample_root_dir` on the oncoanalyser 3.0 output, OA_V3_1 default: no `viridian/` there, so `viridian` null | 6.45 s, 687 MB | 0.27 s, 41 MB | identical, 52,689 bytes |
| hg008_root_format | root with `-pipeline_format OA_V2_0` (virusinterpreter/) and plot directories | 7.92 s, 751 MB | 0.27 s, 42 MB | identical, 53,876 bytes |
| hg008_sample_data_dir | `-sample_data_dir` (the Virus Interpreter file in it is read) | 9.07 s, 749 MB | 0.28 s, 41 MB | identical, 53,817 bytes |
| hg008_bad_format | invalid `-pipeline_format` | 0.34 s, 90 MB | 0.00 s, 2 MB | both exit 1 |
| hg008_no_plots | a PURPLE plot directory without the plots | 0.89 s, 310 MB | 0.21 s, 33 MB | both exit 1 |
| hg008_bad_type | invalid `-experiment_type` | 0.36 s, 93 MB | 0.00 s, 2 MB | both exit 1 |
| hg008_no_tumor | no `-tumor` | 0.32 s, 92 MB | 0.00 s, 2 MB | both exit 1 |
| hg008_wrong_reference | a reference id not in the VCFs | 0.42 s, 150 MB | 0.01 s, 3 MB | both exit 1 |
| t0 | the private pair (renamed copies), the pipeline's command | 7.58 s, 654 MB | 0.17 s, 42 MB | identical, 36,548 bytes |
| t0 :json | the same, harness | 0.88 s, 268 MB | 0.17 s, 42 MB | identical, 36,548 bytes |
| t0_tumor_only | the private pair, tumor-only | 6.13 s, 528 MB | 0.12 s, 33 MB | identical, 20,231 bytes |
| testrun_wgts | the Java tests' run at ca04e6c39d (`virusinterprtr/`, header only) | 6.55 s, 944 MB | 0.55 s, 39 MB | identical, 39,674 bytes |
| testrun_wgts :json | the same, harness | 1.20 s, 615 MB | 0.56 s, 39 MB | identical, 39,674 bytes |
| testrun_panel | its targeted tumor-only run | 5.16 s, 890 MB | 0.49 s, 30 MB | identical, 24,869 bytes |
| testrun_rna_tumor_id | ISOFOX without `-rna_sample` | 1.90 s, 688 MB | 0.55 s, 36 MB | both exit 1 (OR10) |
| testrun_rna_tumor_id :json | the same, harness | 1.15 s, 604 MB | 0.55 s, 39 MB | identical, 22,610 bytes |
| synthetic | the synthetic case with RNA (4 Virus Interpreter viruses, one blacklisted) | 2.71 s, 436 MB | 0.27 s, 41 MB | both exit 1 (OR11) |
| synthetic :json | the same, harness (OR12: the blacklisted virus is listed) | 1.02 s, 341 MB | 0.27 s, 40 MB | identical, 68,703 bytes |
| synthetic_tumor_only :json | tumor-only, harness | 0.93 s, 296 MB | 0.22 s, 32 MB | identical, 48,336 bytes |
| synthetic_no_rna | without ISOFOX | 11.49 s, 894 MB | 0.27 s, 39 MB | identical, 66,011 bytes |
| hg008_format_file | root with the older `-pipeline_format_file` (`virusInterpreterDir`, now ignored) | 1.00 s, 344 MB | 0.28 s, 43 MB | both exit 1 |
| hg008_viridian | NEW: the pipeline's command with `-viridian_dir viridian`: HG008's Virus Interpreter file converted (header only) | 9.40 s, 754 MB | 0.29 s, 39 MB | identical, 53,903 bytes |
| hg008_viridian_dir_vi | NEW: `-viridian_dir virusinterpreter` (a Virus Interpreter file through the new option) | 9.23 s, 758 MB | 0.28 s, 39 MB | identical, 53,903 bytes |
| hg008_no_virus | NEW: the pipeline's command without a virus directory (`viridian` null) | 9.54 s, 768 MB | 0.28 s, 42 MB | identical, 53,852 bytes |
| synthetic_viridian | NEW: synthetic with RNA, `-viridian_dir viridian` (its 4 viruses converted to a Viridian report) | 2.72 s, 436 MB | 0.27 s, 44 MB | both exit 1 (OR11) |
| synthetic_viridian :json | NEW: the same, harness; equals synthetic :json byte for byte | 1.03 s, 344 MB | 0.27 s, 41 MB | identical, 68,703 bytes |
| synthetic_viridian_no_rna | NEW: the same without ISOFOX; equals synthetic_no_rna | 11.31 s, 898 MB | 0.27 s, 42 MB | identical, 66,011 bytes |
| synthetic_viridian_extra | NEW: a hand-made Viridian report: 7 viruses, every clonality status and virus type, other column order, NULL/null, TRUE/False, +3, NaN, Infinity, 1e-7, 12345678.9 | 10.97 s, 892 MB | 0.26 s, 35 MB | identical, 67,588 bytes |
| synthetic_viridian_extra_tumor_only | NEW: the same, tumor-only | 9.91 s, 751 MB | 0.22 s, 33 MB | identical, 47,221 bytes |
| synthetic_virus_dir_viridian | NEW: `-virus_dir` on a directory with only the Viridian report | 11.34 s, 899 MB | 0.26 s, 37 MB | identical, 67,588 bytes |
| synthetic_both_options | NEW: `-virus_dir virusinterpreter -viridian_dir viridian_extra` (`-virus_dir` wins) | 11.43 s, 884 MB | 0.27 s, 42 MB | identical, 66,011 bytes |
| synthetic_both_files | NEW: `-viridian_dir` on a directory with both files (Virus Interpreter's read) | 11.55 s, 892 MB | 0.27 s, 42 MB | identical, 66,011 bytes |
| synthetic_viridian_bad | NEW: a report with the interpretation `HHV8` (IllegalStateException) | 0.94 s, 336 MB | 0.26 s, 38 MB | both exit 1 |
| synthetic_viridian_empty | NEW: `-viridian_dir` on a directory without a report (FileNotFoundException) | 0.93 s, 341 MB | 0.26 s, 43 MB | both exit 1 |
| synthetic_viridian_missing | NEW: `-viridian_dir` that does not exist (invalid path for config) | 0.35 s, 91 MB | 0.00 s, 2 MB | both exit 1 |
| hg008_root_viridian | NEW: a root with `viridian/` (the hand-made report), OA_V3_1 default | 6.73 s, 692 MB | 0.28 s, 37 MB | identical, 56,594 bytes |
| hg008_root_viridian_oa_v3_0 | NEW: the same root with `-pipeline_format OA_V3_0` (virusinterpreter/) | 6.49 s, 684 MB | 0.28 s, 44 MB | identical, 52,740 bytes |
| hg008_format_file_v2 | NEW: the same root with a format file naming `viridianDir` | 6.53 s, 684 MB | 0.28 s, 43 MB | identical, 56,552 bytes |
| testrun_viridian | NEW: the Java tests' run at 883ab11bba, `-viridian_dir virus` (header-only report) | 6.39 s, 944 MB | 0.55 s, 38 MB | identical, 39,674 bytes |
| testrun_virus_dir_viridian | NEW: the same with `-virus_dir virus` | 6.44 s, 946 MB | 0.56 s, 37 MB | identical, 39,674 bytes |

Cross-checks: the Viridian route gives the Virus Interpreter route's JSON byte for byte (`synthetic_viridian :json`
= `synthetic :json`, `synthetic_viridian_no_rna` = `synthetic_no_rna`, `hg008_viridian` = `hg008`), and the 38
identical cases have identical `plot/` folders.

The pipeline's own HG008 JSON (Hartwig's ORANGE image 5.1-beta.1, `orange/HG008-T.orange.json`) now equals both the
jar's and the port's (`-java_compat`, the pipeline's command without `-sampling_date`, run on 2026-11-27, the
pipeline's sampling date) but for the order of the five LINX driver plots (OR3): the image's list is in file name
order (cluster 132, 144, 150, 170, 4), this Mac's directory order is 144, 4, 170, 150, 132. With the entries put in
the image's order the files are byte-identical (53,903 bytes). In the default mode the port sorts them (OR3), which
here is the image's order; the default JSON then differs from the image only in the corrections OR2 (15 plot paths)
and OR5 (`simpleDups32To200B` 5 to 1).

**Default mode** (every case without `-java_compat`; rerun for 9 cases against the new jar, the `viridian` section
equal to Java's in each) against Java differs only where a correction applies: the plot paths (OR2) in
every case; the LINX plot order (OR3) where the directory order is not sorted (HG008); `simpleDups32To200B` (OR5:
HG008 5 to 1, the private pair and the test run too); the DOIDs (OR1); `rnaId` (OR6); in the synthetic case the
other transcripts' exon, codon and coding effect (OR7), the fusions' junction copy numbers (OR8) and the SAGE plot of
a variant whose name is a prefix of another's (OR4). OR9 changes none of these samples (tested).

**Speed.** Java's run is mostly the PDF (PDFBox drawing the large PNGs): without it, 0.9-1.2 s. The port takes
0.12-0.56 s and 30-44 MB (against the new jar; the Viridian report adds nothing measurable), 3-5 times faster than Java's JSON part and 20-40 times faster than the whole jar, in a
tenth of the memory or less. Of the port's time, 85-93% is reading PURPLE's somatic VCF (HG008: 17,254 records,
0.22 s of 0.28 s; decompression is 0.02 s, the rest the parse of each PASS record's INFO and genotypes); LINX's files
take 0.03 s, everything else (interpretation, plot copies, JSON) about 0.01 s. Filtered records are skipped before
parsing in the PASS-only somatic read (as sigs-rs does); a malformed filtered record would fail in Java and not here.
These are the times of the port before it wrote the PDF (now the same with `-no_pdf`); with the PDF, see section
"PDF".

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
  JSON is written either. Corrected: "NA" as the HR deficiency score. Java mode fails at the same point (exit 1, no
  file written).
- OR11: with RNA, a run where one reported DNA fusion has RNA support and another has none fails: the PDF's fusion
  table adds the RNA column once any fusion has support and then reads `rnaSupport()` of every fusion
  (DnaFusionTable.build). No JSON is written. Corrected: "NA" in the RNA Frags column for a fusion without RNA support.
  Java mode fails the same way.
- OR12 (upstream f5f660cc9f; not corrected, Java's behaviour in both modes): Virus Interpreter's blacklisted viruses
  are no longer left out. ORANGE 5.1 at ca04e6c39d dropped them from both lists (VirusInterpreter
  .filterBlacklistedViruses, which also failed on a file without the `blacklisted` column); the new conversion to
  Viridian's reports (VirusInterpreterDataLoader.convertToViridianReports) takes every annotated virus, so a
  blacklisted virus now appears in `allViruses` (Virus Interpreter never reports one: it throws for a virus both
  reported and blacklisted, so `reportableViruses` is unchanged on real files). Seen in the synthetic case, whose
  third virus is marked blacklisted (and, impossibly for Virus Interpreter, reported): it is in both lists now.
  The commit says Virus Interpreter's file loads "as before", which suggests the drop was not intended, but Viridian
  itself has no blacklist and lists every virus group it profiles, so the intent is unclear: kept as Java.
- OR13: in the germline variants table a variant without a gnomAD frequency shows "  nu": the column is
  `String.format("%4.2e", gnomadFrequency())`, and Formatter prints a null argument as "null" cut to the precision (2)
  and padded to the width (4) (FormatCheck harness: `String.format("%4.2e", (Double) null)` is "  nu"). No case has
  such a variant (the germline variants of the cases have frequencies). Corrected: "NA".
- OR14: the scaled tables (variants, gains and deletions, fusions, disruptions) size their columns to their widest text
  (TableCommon.scaleTableWidths, aiming at 450 pt) and stretch them to the table's width (530 pt), but leave out the
  cells' padding (2 + 2 pt). When the text is wider than about 450 pt, the columns shrink below their text plus
  padding, and boxable breaks values inside their cells: in the synthetic cases (Java's own PDFs) the germline allele
  frequencies are drawn as "1." over "00" and "0." over "00", and the fusion table's headers as "JC" over "N" and
  "DRIVE" over "R". Corrected, only for a table where Java breaks a value (a column that may not wrap is narrower than
  its text plus padding): widths of text plus padding, the space left shared as the text widths, the wrappable
  columns (HGVS, a header with WrapAtPercent) giving up width first. HG008's and the private pair's tables are
  unchanged.
- OR15: a title is left alone at the end of a page. A table's title is an ordinary first row before the header row,
  so boxable's rule that keeps the header rows with the first data row (getMinimumHeight, which here counts the header
  row twice) moves the header and rows to the next page and leaves the title (private pair, page 2: the Signature
  Allocations title, with its table on page 3); likewise the "Structural Driver Plots (n)" title, when the plots move
  to a new page because the first row does not fit (synthetic page 4, private pair tumor-only page 2). Corrected: a
  table starts on a new page when its title, header row and first data row (or its title and one row) do not fit
  together; the plots' title moves with the first row of plots.
- OR16: the pages in the middle of a table over three or more pages have no header (ORANGE title, logo) and no side
  panel (the sample id): DocumentContext.addTable runs the page events only for the page the table ends on, and
  addTableNoAdvance (the front page's two tables side by side) for none. No case has such a table; the Java harness
  LayoutStress (a table of 150 rows, a row taller than a page) gives five such pages of nine. Corrected: every page a
  table adds gets the header and side panel.

Not corrected (a data model change, or the intent is unclear), as Java:
- A germline aberration TERTRASOMY_9 in PURPLE's QC (hmf-common has it, the data model's PurpleGermlineAberration
  does not), and germline amp/del statuses LIKELY_DIPLOID, EXCLUDED and CENTROMETIC (not in PurpleGermlineStatus), end
  the run with "No enum constant"; so does an Isofox splice junction of type CANONICAL. The port fails the same way.
- Without a pipeline version file `pipelineVersion` is "OA_V3_1", the name of the default pipeline directory layout.
- The SAGE plot directory's default (with `-pipeline_sample_root_dir`) is under LINX's somatic directory.
- A tool missing from a `-pipeline_format_file` gets the directory "" and so the root directory itself: ORANGE then
  reads e.g. `<root>/<tumor>.cuppa.vis_data.tsv` and fails, or warns "invalid RNA sample ID" for ISOFOX. Since
  f5f660cc9f the virus directory's key is `viridianDir`; a file of the older format names `virusInterpreterDir`,
  which is ignored, so the root directory is taken and the run fails on `<root>/<tumor>.viridian.report.tsv`
  (case `hg008_format_file`).
- Viruses from a Virus Interpreter file (f5f660cc9f's conversion): the virus group is the virus name, every QC
  status other than NO_ABNORMALITIES (LOW_VIRAL_COVERAGE, EXCESSIVE_VIRAL_COVERAGE, ASSEMBLY_DOWNSAMPLED,
  CHILD_TAXID_REFERENCE, UNCLEAR_TAXID_ASSIGNMENT) becomes AMBIGUOUS_VIRUS (Viridian's "virus doesn't match between
  read profiling and integrations"), the clonality is NOT_PROFILED (Viridian's "not present in sample reads") with
  no copies per tumor cell, and the expected clonal coverage and driver likelihood are dropped. The coverage becomes a
  fraction as `percentageCovered * 0.01`; that differs from `/ 100` in the last digit for about 14% of values, but
  neither form prints the decimal in every case (0.47% gives 0.0047 one way and 0.004699999999999999 the other), so
  it is not a fault.
- The default layout (OA_V3_1) now takes the virus directory `viridian/`: a `-pipeline_sample_root_dir` run on an
  oncoanalyser 3.0 output, which has `virusinterpreter/`, gets no `viridian` section (null), with nothing logged above
  TRACE, unless `-pipeline_format OA_V3_0` is given (cases `hg008_root`, `hg008_root_viridian_oa_v3_0`).
- The file read is chosen by what the directory holds, not by the option: a Virus Interpreter file is read whenever
  the virus directory has one, also for `-viridian_dir`; `-virus_dir` reads a Viridian report when there is no Virus
  Interpreter file; with both options `-virus_dir` is used and `-viridian_dir` ignored, silently.
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
- PDF: when the front page's two side-by-side tables do not fit on the page, the left one continues on a new page and
  the right one is then drawn on that page, from the height it had on the first (LayoutStress only; the front page
  has room for them in every case). Kept.
- PDF: the counts in the titles are not the rows below them: "Small Variants (n)" counts the reported transcript
  impacts while the table lists every impact of each variant (canonical and other, reported or not); "Signature
  Allocations (n)" counts all allocations while the table lists those of at least 8% and 8 variants; "Disruptions (n)"
  counts structural variants (a row can show both breakends of one). The intent is unclear; kept.
- PDF: the immune escape table is commented out in ImmunologyChapter (so not drawn); HLA class II alleles only for
  whole genomes; `purpleVariantCopyNumberPlot` (somatic.png) is copied but not drawn.
- PDF: a character Nimbus Sans L has no glyph for (a tab, most non-Latin letters) ends the run in PDFBox
  (IllegalArgumentException "No glyph for U+...") and in the port. Text that boxable reads as markup (`<b>`, `<i>`,
  `<p>`, `<br>`, lists, headings) would be formatted by Java (in Helvetica); the port draws it as text and logs a
  warning (no data has it). A null `biallelicProbability` or, tumor-only, `somaticLikelihood` would end Java's run
  (unboxing); PURPLE always gives both, and the port fails the same way.
- Other failures inside Java's PDF writer than OR10 and OR11 are not reproduced (none was met).

## Tests

- 49 tests in `orange-rs` (`cargo test --release`): the Java unit tests of the ported classes, and tests of the port.
  For the PDF: PlotPathResolverTest, TableDataTest (scaleTableWidths), SignatureAllocationTableTest, MiscUtilsTest
  (zeroPrefixed), ExpressionsTest, QcStatusInterpretationTest and ReportWriterTest (targeted, tumor-only and
  tumor/normal reports from the Java tests' run, in memory); the fonts' metrics and every character width checked
  against PDFBox's (FontMetricsDump harness), the number formats against Java's (FormatCheck harness: DecimalFormat
  "#'%'", "0'%'", "0.0", "0.00", String.format "%.1f", "%.0f", "%4.2e" with null), and `layout_stress_pages`, the
  synthetic layout the Java harness LayoutStress draws (written to `LAYOUT_STRESS_PDF` for the comparison).
  For f5f660cc9f: ViridianInterpreterTest (its two tests), PipelineToolDirectoriesTest.canLoadFromCompleteFile with
  the layouts' virus directories, the VirusReport JSON as `minimally.populated.orange.json` writes it, and the
  loading of either tool's file (both files, neither, no directory). VirusInterpreterTest and OrangeConversionTest's
  virus test were removed upstream, and so here.
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
  viruses, Viridian's virus reports (VirusReport, `read_virus_reports` as VirusReportFile.read,
  `virus_report_filename`, the ReportQcStatus and ClonalityStatus names) and `convert_to_viridian_reports`
  (VirusInterpreterDataLoader.convertToViridianReports), CUPPA's vis_data, PipelineVersionFile.
- `delim_file`: `get_bool` (DelimFileReader.Row.getBoolean: true / false in any case, anything else an error),
  `get_double_or_null`, and `parse_delim_bool`.
- `purple_files`: `read_purity_context_with_qc` (the PurplePurity reader qsee-rs added is shared), ChrArmCopyNumbersFile,
  GermlineAmpDel, `human_chromosome` (HumanChromosome.fromString), the germline driver catalog file name.
- `hla`: LilacQcData (`read_lilac_qc`, `lilac_qc_filename`).
- `rna`: RnaStatisticFile, GeneExpressionFile, RnaFusionFile and NovelSpliceJunctionFile readers.
- `lib.rs`: the new modules and their line in the crate documentation.

## PDF

### How it is drawn

`src/report/` follows `orange/report` (PDFBox 3.0.6 with boxable 1.8.2; about 5,500 lines) class by class:
DocumentContext (cursor, spacing, page breaks, images scaled to fit, never above their size), PageEventHandler with
Header (logo, "ORANGE Report" in six orange shades, "(Research Use Only)" with `-add_disclaimer`), SidePanel (the
sample id, `-display_sample` replacing the tumor id, shrunk by 0.1 pt until 120 pt wide), Footer (page numbers in a
second pass, the disclaimer), the chapters and the tables, from the record and the copied plots only (as Java).

The PDF itself is written by `hmf_plot::pdf` (hmf-plot-rs, feature "pdf", pure Rust): the tables are laid out by its
port of boxable (`pdf::table`), with the metrics PDFBox computes for ORANGE's fonts (Nimbus Sans L Regular and Bold,
embedded as subsets: widths rounded to 1000 units, ascent, descent and x-height from the font tables, boxable's line
height (ascent - descent) / 1000), in the same float arithmetic: the same column widths, line breaks (wrap points after
space , . - @ : ;, words wider than their cell broken), row heights, page breaks with the header row repeated, rows
split across pages. The content streams use PDFBox's operators and number format (5 decimals) and, in tables,
boxable's text handling (one text object per run of rows, each token placed at an offset from the one before), so a
viewer computes the same positions as from Java's file. Images are the PNGs as PDFBox embeds them (a soft mask for
any transparency; 8-bit RGB data passed through without recompression).

### How it was checked

Each case ran the port in Java mode and compared its PDF with Java's (the jar of 883ab11bba): page count and
sizes, `pdftotext -layout` page by page (and the words without spacing), every word's box (`pdftotext -bbox`),
the images (`pdfimages -list`: page, pixel size, colour space, soft masks), the JSON, and the pages rendered at
100 dpi (`pdftoppm`), shown side by side and their pixels compared (ImageMagick, fuzz 2%).
Three Java harnesses gave references beyond the cases: PDFBox's font metrics (FontMetricsDump), Java's number
formats (FormatCheck), and a layout stress test drawn with ORANGE's own table classes (LayoutStress).

Java mode, the 29 cases where Java writes a PDF (time and peak RSS of the whole run, `/usr/bin/time -l`, arm64 Mac):

| Case | Pages | Words (all at the same place) | Images | Java | Rust |
|---|---|---|---|---|---|
| hg008 | 9 | 1,078 | 34 | 9.14 s, 778 MB | 0.65 s, 131 MB |
| hg008_doids | 9 | 1,077 | 34 | 9.55 s, 742 MB | 0.66 s, 138 MB |
| hg008_format_file_v2 | 8 | 813 | 27 | 6.42 s, 678 MB | 0.64 s, 135 MB |
| hg008_location | 9 | 1,077 | 34 | 9.38 s, 751 MB | 0.65 s, 130 MB |
| hg008_minimal (no disclaimer) | 5 | 524 | 19 | 4.12 s, 541 MB | 0.57 s, 110 MB |
| hg008_no_virus | 9 | 1,072 | 34 | 9.27 s, 741 MB | 0.66 s, 131 MB |
| hg008_options (`-display_sample`, `-output_id`, SBX ...) | 9 | 1,084 | 34 | 9.19 s, 753 MB | 0.65 s, 134 MB |
| hg008_panel_tumor_only | 7 | 827 | 29 | 7.70 s, 696 MB | 0.60 s, 128 MB |
| hg008_root | 8 | 768 | 27 | 6.45 s, 675 MB | 0.65 s, 139 MB |
| hg008_root_format | 8 | 792 | 31 | 8.20 s, 738 MB | 0.65 s, 139 MB |
| hg008_root_viridian | 8 | 813 | 27 | 6.37 s, 676 MB | 0.65 s, 135 MB |
| hg008_root_viridian_oa_v3_0 | 8 | 774 | 27 | 6.49 s, 691 MB | 0.64 s, 133 MB |
| hg008_sample_data_dir | 9 | 799 | 34 | 9.17 s, 757 MB | 0.66 s, 143 MB |
| hg008_tumor_only | 7 | 890 | 29 | 7.73 s, 689 MB | 0.61 s, 123 MB |
| hg008_viridian | 9 | 1,078 | 34 | 9.38 s, 755 MB | 0.68 s, 139 MB |
| hg008_viridian_dir_vi | 9 | 1,078 | 34 | 9.51 s, 765 MB | 0.66 s, 134 MB |
| synthetic_both_files | 12 | 1,440 | 44 | 11.29 s, 897 MB | 0.65 s, 130 MB |
| synthetic_both_options | 12 | 1,440 | 44 | 11.25 s, 896 MB | 0.65 s, 131 MB |
| synthetic_no_rna | 12 | 1,440 | 44 | 11.59 s, 878 MB | 0.67 s, 119 MB |
| synthetic_viridian_extra | 12 | 1,452 | 44 | 11.47 s, 892 MB | 0.66 s, 141 MB |
| synthetic_viridian_extra_tumor_only | 8 | 1,039 | 35 | 10.05 s, 754 MB | 0.59 s, 134 MB |
| synthetic_viridian_no_rna | 12 | 1,440 | 44 | 11.55 s, 888 MB | 0.65 s, 141 MB |
| synthetic_virus_dir_viridian | 12 | 1,452 | 44 | 11.12 s, 876 MB | 0.66 s, 129 MB |
| t0 (the private pair) | 8 | 847 | 29 | 7.39 s, 655 MB | 0.54 s, 138 MB |
| t0_tumor_only | 6 | 654 | 24 | 6.02 s, 529 MB | 0.49 s, 116 MB |
| testrun_panel | 6 | 532 | 23 | 5.17 s, 899 MB | 0.88 s, 143 MB |
| testrun_viridian | 8 | 637 | 28 | 6.27 s, 896 MB | 0.95 s, 150 MB |
| testrun_virus_dir_viridian | 8 | 637 | 28 | 6.24 s, 913 MB | 0.95 s, 144 MB |
| testrun_wgts | 8 | 637 | 28 | 6.31 s, 913 MB | 0.96 s, 145 MB |

In every one: the same page count and page sizes, `pdftotext -layout` output identical, every word at the same place
(largest difference 0.0000 pt), the same images on the same pages (pixel size, colour space, soft masks), the rendered
pages identical, the JSON byte-identical, the same chapter bookmarks and document title and author. The five cases
where Java fails in its PDF writer (hg008_panel, hg008_minimal_reference, testrun_rna_tumor_id: OR10; synthetic,
synthetic_viridian: OR11) fail in Java mode at the same point, with nothing written. The synthetic layout of the Java
harness LayoutStress (ORANGE's own DocumentContext, Tables and Cells: a scaled table of 150 rows with wrapped and broken
words over four pages, a row taller than a page, a wrapped header, side-by-side tables across a page break, empty
tables to the page end) is equal in the same way (9 pages, 2,960 words).

Not the same: the files' bytes (PDFBox's object layout, compression and subset tags; the port's HG008 file is 10.2 MB
against 9.2 MB, as decoded images are compressed at a faster level), and the PDF version is written as PDFBox's (1.6).

Default mode: the PDFs differ from Java's only by the corrections: the LINX plots in file name order (OR3), the fusion
junction copy numbers (OR8), the configured primary tumor (OR1, test run), a title moved to the next page with its
table or plots (OR15: private pair, synthetic), the synthetic tables' values and headers no longer broken (OR14), and
PDFs where Java fails (OR10: "NA" as HR deficiency score; OR11: "NA" as RNA fragments). The default-mode HG008 PDF
equals the one Hartwig's ORANGE image (5.1-beta.1) wrote in the pipeline (`data/public/hg008/pipeline_full/orange`):
the same text, word positions, images and rendered pages (the image lists the LINX plots in file name order, as OR3
does); in Java mode the order is this Mac's directory order (OR3).

Determinism: two runs give the same bytes, and the Linux arm64 build (`rust:1.97-bookworm`, without CMake) gives the
macOS one's bytes for HG008.

### Speed

The PDF takes 0.36 s of the port's 0.65 s for HG008 (reading the data takes 0.26 s): 0.34 s for the pages, nearly all
of it the images (reading 8 MB of PNGs, and decoding the three with transparency for their soft masks: the header logo
and two 3000 x 3000 Circos plots), 4 ms to make the file, 14 ms to write its 10 MB. Java takes 4-12 s and 530-910 MB
for the whole run, mostly in its PDF writer (without it, about 1 s); the port 0.5-1.0 s and 110-150 MB: 6-18 times
faster, in 13-22% of the memory. Both run on one thread. The port's memory is mostly the plots held for the file
(their data is written at the end) and the decoded images with transparency while their masks are made.

## Open points

- The PDF's fonts and header image are embedded from `orange-rs/resources/` (copied from `orange/src/main/resources`):
  the image's build context needs `!orange-rs/resources` in `rust-docker/Dockerfile.dockerignore`.
- In the pipeline, oncoanalyser runs ORANGE only when PURPLE's `plot/<tumor>.circos.png` exists. Until purple-rs draws
  its plots, a run with purple-rs skips ORANGE, and orange-rs runs only beside the Java PURPLE.
- `small_variant` in hmf-common-rs repeats code that purple-rs has in `src/variant.rs` (VariantImpact, the pathogenic
  summary); purple-rs could use the shared module.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 12 on Google Batch ran ORANGE as the port, in Java mode, on the inputs of rerun 11 (every other step
from the cache; PURPLE in Java, which draws the plots ORANGE needs). The JSON equaled the one Hartwig's ORANGE image
(5.1-beta.1) wrote from the same inputs but for one key: the image names the virus section `viridian`, while the
port, then following ca04e6c39d, named it `virusInterpreter`. The port now follows f5f660cc9f; run locally on the
same inputs with the pipeline's command it gives the image's JSON but for the LINX plot order (above). The LINX
driver plots came in the same order on Batch (in Java mode both take Linux's directory order). ORANGE took 5.3 s and
42 MB, against Java's 18 s and 419 MB.

**The PDF (HG008 rerun 15):** orange-rs with its PDF ran in Java mode on rerun 11's inputs (Java PURPLE and its R and
Circos plots). Against the PDF Java ORANGE wrote from the same inputs: 8 pages of the same sizes, the same
`pdftotext -layout` text, every one of 1,041 words at the same position, the same 30 embedded images, no rendered
pixel different at 100 dpi; the JSON byte-identical. 5 s and 43 MB, against Java's 18 s and 419 MB.
