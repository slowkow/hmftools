# orange-rs

A Rust port of hmftools ORANGE 5.1 (`../orange`, with the data model `../orange-datamodel`): the report that gathers
every tool's results for a tumor sample. It takes ORANGE's command line, loads and interprets PURPLE, LINX, LILAC,
CHORD, CUPPA, PEACH, SIGS, Virus Interpreter and ISOFOX results as ORANGE does (driver and variant selection,
reportable events, germline handling, immune escape), copies the plots, and writes `<tumor>.orange.json`. With
`-java_compat` the JSON is byte-identical to the Java jar's; by default the faults found in the Java code are
corrected (findings OR1-OR11 in REPORT.md).

The PDF report is not written (PDFBox and boxable in Java; out of scope for now). REPORT.md lists each section and
table the PDF draws and the data and plots it uses.

## Build and run

```
CARGO_BUILD_JOBS=2 cargo build --release          # target/release/orange-rs
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
  `-cuppa_dir`, `-peach_dir`, `-sigs_dir`, `-virus_dir`, `-qsee_dir`, `-isofox_dir`, `-purple_rna_append_vcf`; or
  `-pipeline_sample_root_dir` with `-pipeline_format` (OA_V3_1 by default, OA_V3_0, OA_V2_3, OA_V2_2, OA_V2_0,
  PIP5_V6_0, DB_V6_0) or `-pipeline_format_file`; or `-sample_data_dir` (everything in one directory).
- **Report:** `-pipeline_version_file`, `-sampling_date` (yyMMdd; today by default), `-primary_tumor_doids` (with
  `-doid_json`), `-primary_tumor_location`, `-panel_name`, `-display_sample`, `-add_disclaimer` (the last two are
  PDF-only).
- **Output:** `-output_dir`, `-output_id` (`<tumor>.orange.<id>.json`), `-log_level`, `-log_debug`.
- **Port:** `-no_pdf` (accepted; the port never writes the PDF), `-java_compat` (Java's behaviour for every
  finding) and `-java_compat_except OR1,...` (all but those listed); the environment variables
  `ORANGE_JAVA_COMPAT=1` and `ORANGE_JAVA_COMPAT_EXCEPT` do the same.

## Outputs

In `-output_dir`:
- `<tumor>.orange.json`: the report, as Gson writes it (pretty-printed, HTML-safe escapes, nulls written);
- `plot/`: the PURPLE plots, the LINX driver plots, the SAGE variant pages, the QSee and CUPPA plots, copied.

## How it works

`config.rs` resolves the options and directories (OrangeConfig, PathResolver, the pipeline layouts); `purple.rs`,
`linx.rs`, `isofox.rs` and `others.rs` load each tool's files (through hmf-common-rs: `purple_files`,
`small_variant`, `linx_files`, `rna`, `hla`, `tool_results`, `cytobands`) and interpret them as OrangeAlgo and its
interpreters do; `model.rs` is the data model, each type writing its attributes in its interface's order as the
Immutables Gson adapters do, through `json.rs` (Gson's JsonWriter). Java's HashSet order is reproduced where it reaches
the JSON (the local phase set ids), EnumSet order for enum sets, and directory listing order for the plots (OR3).

## Tests

`CARGO_BUILD_JOBS=2 cargo test --release` runs the Java unit tests of the ported classes (GermlineGainDeletionFactory,
PurpleInterpreter, CodingEffectDeterminer, DriverInterpretation, LinxBreakendInterpreter, ImmuneEscapeInterpreter,
Isofox interpretation and conversion, VirusInterpreter, CuppaDataFactory with the Java tests' CUPPA files,
FileBasedPlotManager, PathUtil) and tests of the port (JSON layout, findings in both modes, dates, plot lookup).
The comparisons with Java are in REPORT.md.
