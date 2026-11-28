# orange-rs

A Rust port of hmftools ORANGE 5.1 (`../orange`, with the data model `../orange-datamodel`; upstream 883ab11bba of
2026-11-27, which has f5f660cc9f's Viridian input): the report that gathers every tool's results for a tumor sample.
It takes ORANGE's command line, loads and interprets PURPLE, LINX, LILAC, CHORD, CUPPA, PEACH, SIGS, Viridian (or
Virus Interpreter) and ISOFOX results as ORANGE does (driver and variant selection, reportable events, germline
handling, immune escape), copies the plots, and writes the PDF report `<tumor>.orange.pdf` and `<tumor>.orange.json`.
With `-java_compat` the JSON is byte-identical to the Java jar's and the PDF has Java's pages: the same text at the
same positions (to the last digit a viewer reads), the same tables, page breaks, plots and page images; by default
the faults found in the Java code are corrected (findings OR1-OR16 in REPORT.md; OR12 is kept as Java).

The PDF is drawn natively with the PDF backend of `hmf-plot-rs` (pure Rust): Java's PDFBox and boxable layout is
followed step by step, with ORANGE's fonts (Nimbus Sans L) and header image, which are embedded in the binary
(`resources/`).

## Build and run

```
CARGO_BUILD_JOBS=2 cargo build --release          # target/release/orange-rs (pure Rust apart from hmf-common-rs's htslib; no CMake)
orange-rs -experiment_type WGS -tumor TUMOR -reference NORMAL -purple_dir purple -purple_plot_dir purple/plot/ \
    -linx_dir linx_somatic -linx_germline_dir linx_germline -linx_plot_dir reportable/ -sage_plot_dir sage_vis \
    -lilac_dir lilac -chord_dir chord -cuppa_dir cuppa -peach_dir peach -sigs_dir sigs -virus_dir virusinterpreter \
    -qsee_dir qsee -pipeline_version_file pipeline_version.txt -ref_genome_version 38 -output_dir output/
```

This is oncoanalyser's ORANGE task (its `-add_disclaimer`, `-sequencing_type`, `-sage_dir`, `-sage_germline_dir` and
`-doid_json` are accepted as well). JVM options (`-Xmx...`, `-Xms...`, `-XX:...`, `-D...=...`) are accepted and
ignored, so the bioconda wrapper's command line works, and `orange-rs -version` prints `Orange version 5.1`.

## Options

As Java's (`orange-rs -help` lists them):
- **Run:** `-experiment_type` WGS or PANEL (required), `-tumor` (required), `-reference` (tumor/normal; without it
  tumor-only), `-rna_sample`, `-ref_genome_version` (37 by default, or 38), `-sequencing_type` (ILLUMINA, ULTIMA, SBX).
- **Inputs:** `-purple_dir`, `-purple_plot_dir`, `-linx_dir` (required, or resolved from a root directory),
  `-linx_germline_dir` (required with a reference), `-linx_plot_dir`, `-sage_plot_dir`, `-lilac_dir`, `-chord_dir`,
  `-cuppa_dir`, `-peach_dir`, `-sigs_dir`, `-viridian_dir`, `-virus_dir`, `-qsee_dir`, `-isofox_dir`,
  `-purple_rna_append_vcf`; or `-pipeline_sample_root_dir` with `-pipeline_format` (OA_V3_1 by default, OA_V3_0,
  OA_V2_3, OA_V2_2, OA_V2_0, PIP5_V6_0, DB_V6_0) or `-pipeline_format_file`; or `-sample_data_dir` (everything in
  one directory).
- **Viruses:** `-virus_dir` (Virus Interpreter's directory, kept for compatibility) is taken as given and wins over
  `-viridian_dir`; without it the virus directory is `-viridian_dir`, else the layout's (`viridian` in OA_V3_1,
  `virusinterpreter` in the older layouts; `viridianDir` in a `-pipeline_format_file`), else `-sample_data_dir`. In
  that directory `<tumor>.virus.annotated.tsv` (Virus Interpreter) is read when it exists, converted to Viridian's
  reports, and otherwise `<tumor>.viridian.report.tsv` (Viridian), which must then exist. No directory, or one that
  does not exist: no `viridian` section (null).
- **Report:** `-pipeline_version_file`, `-sampling_date` (yyMMdd; today by default), `-primary_tumor_doids` (with
  `-doid_json`), `-primary_tumor_location`, `-panel_name`, `-display_sample` (the sample id in the PDF's side panel
  instead of the tumor id), `-add_disclaimer` (the "Research Use Only" title and footer on every PDF page).
- **Output:** `-output_dir`, `-output_id` (`<tumor>.orange.<id>.json`), `-log_level`, `-log_debug`.
- **Port:** `-no_pdf` (only the JSON; Java has no such option, and with `-java_compat` the PDF's failures OR10 and
  OR11 then do not happen), `-java_compat` (Java's behaviour for every finding) and `-java_compat_except OR1,...` (all but those listed); the environment variables
  `ORANGE_JAVA_COMPAT=1` and `ORANGE_JAVA_COMPAT_EXCEPT` do the same.

## Outputs

In `-output_dir`:
- `<tumor>.orange.pdf` (`<tumor>.orange.<id>.pdf` with `-output_id`): the report, A4 pages: Front Page, Somatic
  Findings, Germline Findings (with a reference), Immunology, RNA Findings (with ISOFOX), Tissue of Origin (with a
  reference), Quality Control (with a QSee plot), Purity and Ploidy (landscape for the last three), each with its
  bookmark; written before the JSON, as Java does;
- `<tumor>.orange.json`: the report, as Gson writes it (pretty-printed, HTML-safe escapes, nulls written); the virus
  section is `viridian` (`allViruses` and `reportableViruses`, each a VirusReport: name, virusGroup, qcStatus,
  integrations, interpretation, percentageCovered as a fraction, meanCoverage, reported, clonalityStatus,
  copiesPerTumorCell);
- `plot/`: the PURPLE plots, the LINX driver plots, the SAGE variant pages, the QSee and CUPPA plots, copied.

## How it works

`report/` writes the PDF: `mod.rs` is Java's DocumentContext and page events (header, side panel, footers),
`chapters.rs` the chapters, `tables.rs` the tables and `common.rs` the cells, styles, number formats and column
widths (TableCommon, Cells, Tables); the tables are laid out by `hmf_plot::pdf::table`, a port of boxable's.
`config.rs` resolves the options and directories (OrangeConfig, PathResolver, the pipeline layouts); `purple.rs`,
`linx.rs`, `isofox.rs` and `others.rs` load each tool's files (through hmf-common-rs: `purple_files`,
`small_variant`, `linx_files`, `rna`, `hla`, `tool_results`, `cytobands`) and interpret them as OrangeAlgo and its
interpreters do; `model.rs` is the data model, each type writing its attributes in its interface's order as the
Immutables Gson adapters do, through `json.rs` (Gson's JsonWriter). Java's HashSet order is reproduced where it reaches
the JSON (the local phase set ids), EnumSet order for enum sets, and directory listing order for the plots (OR3).

## Tests

`CARGO_BUILD_JOBS=2 cargo test --release` runs the Java unit tests of the ported classes (GermlineGainDeletionFactory,
PurpleInterpreter, CodingEffectDeterminer, DriverInterpretation, LinxBreakendInterpreter, ImmuneEscapeInterpreter,
Isofox interpretation and conversion, ViridianInterpreter, CuppaDataFactory with the Java tests' CUPPA files,
FileBasedPlotManager, PathUtil, PipelineToolDirectories' complete format file, and of the report: PlotPathResolver,
TableData's width scaling, SignatureAllocationTable's sort, zeroPrefixed, Expressions, QcStatusInterpretation, and
ReportWriterTest's reports from the test run) and tests of the port (the virus loading from either tool's file, JSON
layout, findings in both modes, dates, plot lookup; the fonts' metrics against PDFBox's, the number formats against
Java's, and a synthetic layout that a Java harness draws the same way).
The comparisons with Java are in REPORT.md.

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20), with upstream's Viridian loading (`f5f660cc9f`) taken while porting: ORANGE as at `883ab11bba` (2026-11-27). `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
