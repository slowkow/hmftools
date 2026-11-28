# qsee-rs

A Rust port of hmftools qsee 1.0 (`../qsee`), the QSEE_METRICS step of oncoanalyser: QC features collected
from the PURPLE, BamTools, COBALT, ESVEE and REDUX outputs of a tumor (and its normal), QC statuses against
thresholds, the limit of detection of variants per exon of the reportable driver genes, and a comparison
with a cohort. It takes the same options and writes the same files as the Java `QseeApplication`; with
`-java_compat` they are byte-identical to Java's (gzip bytes included). Faults found in the Java code are
corrected by default (findings QS1-QS16 in [REPORT.md](REPORT.md); QS15 and QS16 are kept as R draws them).

## Build and run

```
cargo build --release          # target/release/qsee-rs
qsee-rs -tumor TUMOR -reference NORMAL \
    -redux_tumor_dir redux_tsvs_tumor/ -redux_ref_dir redux_tsvs_normal/ \
    -bam_metrics_tumor_dir bamtools_tumor -bam_metrics_ref_dir bamtools_normal \
    -cobalt_dir cobalt -esvee_dir esvee -purple_dir purple -sequencing_type ILLUMINA \
    -driver_gene_panel DriverGenePanel.38.tsv -cohort_percentiles_file qsee.cohort.percentiles.tsv.gz \
    -allow_missing_input -log_level INFO -output_dir qsee/
```

This is the command line oncoanalyser runs (QSEE). JVM options (`-Xmx...`, `-XX:...`) are accepted and
ignored, so the bioconda wrapper's command line works, and `qsee-rs -version` prints `Qsee version 1.0`.

The jar's other entry points (`java -cp qsee.jar <class>`) are selected by the class name, full or simple,
as the first argument:
- `CohortPercentilesTrainer`: the cohort percentiles file from many samples (`-sample_id_file`), with
  `-num_percentiles` (101), `-percentile_increment` and `-write_cohort_features`;
- `DefaultThresholdsWriter`: `qsee.thresholds.default.tsv`, with `-targeted_mode` and `-no_overridable_only`;
- `GeneLodCalculator`: the gene limit of detection file only;
- `QseePlot`: the report only, from an existing vis data file (`-vis_data_file`).

## Options

As Java's (`qsee-rs -help` lists them):
- **Samples:** `-tumor` and `-reference` (comma-separated lists), or `-sample_id_file` (a TSV with `TumorId`
  and optionally `ReferenceId`). More than one sample of a type is multi-sample mode: the files are named
  `multisample.qsee.*` and the features of each sample are in key order. Without `-reference`, tumor-only.
- **Inputs:** `-sample_data_dir` (all files in one directory), or per tool `-bam_metrics_tumor_dir`,
  `-bam_metrics_ref_dir`, `-redux_tumor_dir`, `-redux_ref_dir`, `-cobalt_dir`, `-esvee_dir`, `-purple_dir`; a
  `*` in a directory is replaced by the sample id. `-allow_missing_input`: a missing input is a warning and its
  features are left out, rather than an error.
- **Resources:** `-driver_gene_panel` (the reportable genes of the missed-variant likelihood and the LoD),
  `-cohort_percentiles_file` (checked for features the cohort lacks; the report draws its percentiles),
  `-threshold_overrides_file`.
- **Other:** `-sequencing_type` (ILLUMINA, SBX, ULTIMA: the base quality bins of the BQR features),
  `-targeted_mode` (the targeted thresholds), `-skip_categories` (comma-separated), `-write_lod_details`
  (the depth counts in the LoD file), `-skip_plot`, `-output_id`, `-threads` (samples in parallel in
  multi-sample mode), `-log_level`, `-log_debug`; for the report `-show_plot_warnings`, `-merge_plots`.
- **Java behaviour:** `-java_compat` gives Java's behaviour for every finding; `-java_compat_except QS1,...`
  for all but those listed. The environment variables `QSEE_JAVA_COMPAT=1` and `QSEE_JAVA_COMPAT_EXCEPT` do
  the same.

## Outputs

- `<tumor>.qsee.vis.data.tsv.gz`: every feature with a value: sample, sample type, source tool, feature
  type and name, the value (`DecimalFormat("0.########")`, `Inf` for infinity) and the plot metadata (group,
  display name, number format, QC status and threshold). This is what the R report draws.
- `<tumor>.qsee.status.tsv.gz`: the features with a QC status (PASS, WARN, FAIL) and the condition.
- `<tumor>.qsee.gene_lod.tsv.gz`: per exon of the reportable genes, the sample and tumor VAF at which a
  variant is detected 95% of the time, and PASS / WARN_HET_UNRELIABLE / FAIL.
- `<tumor>.qsee.vis.report.pdf` and (for a single patient) `.png`: the report that Java has R draw with
  `plot_qc.R`, drawn natively (`src/plot/`, with hmf-plot-rs), unless `-skip_plot`. With `-merge_plots` in
  multi-sample mode the reports' pages go into `multisample.qsee.vis.report.pdf` and the per-sample PDFs are
  deleted, as Java does with PDFBox.

The gzip files are written as java.util.zip writes them (zlib level 6, OS 255), so their bytes are Java's.

## The report

The port draws what the repository's `plot_qc.R` (copied unchanged into `resources/` as the reference, not
run) draws with R 4.5, ggplot2 4.0.3 and patchwork 1.3.2, from the same inputs as QseePlotTask: tumor id,
reference id or none, the vis data file, the cohort percentiles file or none, the PDF path, single patient.
R is not needed. The layout equals R's (every element within 0.01 px of R's viewports, in the PNG and the PDF);
text in the PNG is DejaVu Sans as in the qsee image but rasterised by hmf-plot-rs (whole-pixel glyph positions),
so pixels at glyph edges, at panel edges (whole-pixel clipping) and of anti-aliased lines differ from Cairo's.
The PDF has R's page, text and word positions (Helvetica with R's metrics and kerning) and draws as R's at
100 dpi to within a few pixels. The legend is sized with PDF metrics in the PNG too, as the script builds it
while R's `pdf(NULL)` device is open. The comparisons are in [REPORT.md](REPORT.md), "The report, drawn
natively". In the image, the `qsee` wrapper must pass `-skip_plot` on only when the caller gives it (the old
wrapper added it whenever Rscript was missing).

What the report draws, and from which data (all from the vis data file, with the cohort percentiles file's
5th, 25th, 50th, 75th and 95th percentiles where given):

| Panel | Feature type | From |
|---|---|---|
| A summary table: mapping, copy number, contamination, mutational burden, with WARN/FAIL boxes | SUMMARY_TABLE | BamTools summary, coverage and flag counts; PURPLE purity and QC |
| B coverage distribution | COVERAGE_DISTRIBUTION | BamTools `.coverage.tsv` |
| C fragment length distribution | FRAG_LENGTH_DISTRIBUTION | BamTools `.frag_length.tsv` |
| D duplicate frequency | DUPLICATE_FREQ | REDUX `.duplicate_freq.tsv` (tumor) |
| E GC bias | GC_BIAS | COBALT `.gc.median.tsv` |
| F discordant fragment frequency | DISCORDANT_FRAG_FREQ | ESVEE `.disc_stats.tsv` (tumor) |
| G genes with missed variants (up to 20 above 2%) | MISSED_VARIANT_LIKELIHOOD | BamTools `.gene_coverage.tsv`, PURPLE purity and sex, the driver panel |
| H BQR by original base quality | BQR_PER_ORIG_QUAL | REDUX `.bqr.tsv` |
| I BQR by SNV96 context | BQR_PER_SNV96_CONTEXT | REDUX `.bqr.tsv` (high qualities) |
| J MS indel error rates | MS_INDEL_ERROR_RATES | REDUX `.ms_table.tsv.gz` |
| K MS indel error bias | MS_INDEL_ERROR_BIAS | REDUX `.ms_table.tsv.gz` |

The title gives the overall status (FAIL if any FAIL, else WARN if any WARN) and every WARN and FAIL.

## How it works

Each of the ten categories reads its tool's files, as Java tells a missing input (an IOException: no such
file, not readable, not UTF-8) from an input it cannot parse (which ends the run, `-allow_missing_input` or
not). The readers are hmf-common-rs's (`metrics`, `purple_files`, and the new `redux_files`,
`cobalt_files`, `esvee_files`, `delim_file`). The missed variant likelihood uses commons-math3's binomial
distribution (hmf-common-rs `commons_math`) and the LoD bisects it 20 times; the BQR means are
commons-math3's weighted Mean with its second-pass correction; the MS indel Phred scores use
StrictMath.log10, which is what Java's Math.log10 is on the aarch64 machine of the Java reference.

## Tests

`cargo test --release` runs the Java unit tests of qsee (16 classes, `src/tests/java_tests.rs`) and tests of
the port (`src/tests/port.rs`: number formats, the threshold display, key order, quality bins per sequencing
type, and the corrected findings in both modes; `src/plot/`: R's merge order of the cohort rows, the label boxes
against R's). The report's layout model is tested in hmf-plot-rs (`ggp`, `patchwork`, `pdf_figure`, `raster`).
The comparisons with Java and R are in [REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
