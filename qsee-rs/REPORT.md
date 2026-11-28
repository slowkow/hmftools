# qsee in Rust: record of the port

This report covers how `qsee-rs` was checked against Java qsee 1.0, the faults found in the Java code, and the
timings. `qsee-rs` is the QSEE_METRICS step of oncoanalyser: QC features from the PURPLE, BamTools, COBALT, ESVEE
and REDUX outputs, QC statuses against thresholds, the exon coverage QC file (the gene limit-of-detection file before
upstream bdceb2302a), and the report. The Java reference was this repository's
`qsee/target/qsee-1.0-jar-with-dependencies.jar` (sources following upstream ca04e6c39d of 2026-11-20), on Java 23;
since 2026-11-28 the port follows upstream master 883ab11bba ("Following upstream" below), whose jar and R script are
the reference of `-java_compat`. `resources/plot_qc.R` is that R script, unchanged; since 2026-11-27 the port draws
its report natively (no R). Options and
outputs are in [README.md](README.md).

## What is ported

All of qsee 1.0's main code, every entry point:

| Java | Rust |
|---|---|
| QseeApplication, QseeConfig, QseePlotConfig, TrainConfig | `main.rs`, `app.rs`, `config.rs` |
| common: SampleType, QseeFileCommon, BinnedFrequencies, LinearInterpolation, MultiFieldStringBuilder, SampleIdsLoader | `common.rs` |
| feature: FeatureType, SourceTool, NumberFormat, FeatureKey (its compareTo), Feature, FeatureMetadata, SampleFeatures; category enums (SummaryTableFeature/Group, DiscordantFragGroup, BaseQualBin, RepeatType) | `feature.rs` |
| feature.FeaturePrep, CategoryPrepTask and the ten CategoryPreps (summary table BamMetrics / PURPLE, coverage and fragment length distributions, missed variant likelihood, duplicate frequency, GC bias, discordant fragments, BQR with BaseQualBinner, MS indel errors) | `prep.rs`, `app.rs` |
| status: ComparisonOperator, QcStatusType, QcStatus (with the BigDecimal display), ThresholdKey, QcThreshold, ThresholdRegistry, ThresholdOverridesFile, PurpleQCStatusConverter, QcStatusFile, DefaultThresholdsWriter | `status.rs` |
| coverage (lod before upstream bdceb2302a): LimitOfDetectionCalcs, ExonCoverageQc, ExonCoverageQcFile(Builder), ExonExclusions(Config, FileBuilder), VariantDetectabilityStatus | `coverage.rs` |
| cohort: FeatureMatrix, PercentileTransformer, FeaturePercentiles, CohortPercentiles(File), CohortPercentilesTrainer, CohortFeaturesWriter | `cohort.rs` |
| plot: VisDataFile, QseePlot (with PDFBox's merging), QseePlotTask; `plot_qc.R` | `vis.rs`; `plot/` (native, hmf-plot-rs) |
| hmf-common: DelimFileReader/Writer, createBufferedWriter (gzip), BqrFile, DuplicateFrequency, JitterCountsTableFile, CobaltGcMedianFile, EsveeDiscordantStats, PurplePurity (full), PurityContextFile.readWithQC | hmf-common-rs additions (below) and `prep.rs` |

Not ported: `compare_cohorts.R`, which no Java class runs.

### The plot step

Until 2026-11-27 the port ran the jar's `plot_qc.R` with the 8 arguments QseePlotTask passes (`Rscript`, the
script from a temporary file, as RExecutor); the R step reads only the vis data file and the cohort
percentiles file, and on the two pipeline cases the reports drawn from the port's files equalled Java's (the
PNG byte for byte, the PDF apart from its CreationDate and ModDate). The rows below marked "R report
included" are from then. The port now draws the report itself; see "The report, drawn natively" below.

## Java behaviour reproduced (not faults)

- **Missing inputs:** Java tells an IOException (no file, a directory, not UTF-8: Files.readAllLines'
  MalformedInputException) from other exceptions. The first is a warning with `-allow_missing_input` (the
  category gives no features) or the end of the run; the second ends the run either way. The port reads each
  input first (`prep::read_input`) and parses after, so the same files fail the same way (cases
  `missing_*`, `edge_not_utf8`). The BamMetrics summary table logs its missing files itself and keeps the
  features of the files there are.
- **Where Java stops:** the data files are written in the order vis data, status, gene LoD, so a missing
  exon coverage file without `-allow_missing_input` leaves the first two written and exits 1
  (`missing_no_exon_coverage`); a line Java cannot format leaves a vis data file closed after the lines before
  it (QS6). The port leaves the same files.
- **Numbers:** `DecimalFormat("0.########")` with "Inf" (vis data, LoD, cohort files) and DelimFileWriter's
  `"#.####"` (status values, thresholds), both from hmf-common-rs `fmt::decimal_format` (FloatingDecimal
  digits, HALF_EVEN). Thresholds are shown as BigDecimal.valueOf(d) (times 100 for percentages),
  stripTrailingZeros, toPlainString. The missed variant likelihood uses commons-math3's
  BinomialDistribution.cumulativeProbability (hmf-common-rs `commons_math`); the BQR changes are
  commons-math3 Mean.evaluate(values, weights) with its correction pass; the MS indel Phred scores use
  `StrictMath.log10` (Java's Math.log10 on aarch64 is fdlibm's; on x86-64 HotSpot has an intrinsic that may
  differ in the last bit, which would show only at the 8th decimal in rare cases); totals are Kahan sums
  where Java's DoubleStream.sum is.
- **Orders:** single-sample features in category order (EnumMaps in enum order, LinkedHashMaps in first-seen
  order, BQR records sorted stably by read type, signed ref and alt bytes and context); multi-sample features
  in FeatureKey.compareTo order (ties, as within BQR bins, kept in first-added order); samples sorted by id
  (stable) in multi-sample mode.
- **Edge inputs Java rejects** (exit 1, reproduced): a distribution file with a header only (BinnedFrequencies'
  System.arraycopy of length -1), a BQR quality past 127 ((byte) cast, no bin), an unknown repeat unit, a
  short discordant stats line, an unknown gender, a cohort file with a "null" value, an
  unsupported or misspelt threshold override, an unknown sequencing type, category or log level, `-threads`
  not a number, sample count mismatches, `-percentile_increment 30`.
- **Gzip:** java.util.zip's GZIPOutputStream (header without time stamp, OS 255, zlib level 6); flate2 with
  zlib gives the same bytes (every gzip file compared byte for byte).

## Cases and results

All against the jar, port in Java mode (`QSEE_JAVA_COMPAT=1`); "identical" means the exit statuses agree and
every data file Java wrote equals the port's byte for byte (gzip bytes included), with the same set of
files. Java with `-Xmx8g` (the pipeline's `-Xmx...` is the wrapper's JVM option; the run script drops it from
the program arguments), the R step only in the two pipeline cases (`-skip_plot` elsewhere). One run at a
time on a shared Mac (load about 4); times are wall times, memory the peak resident set (`/usr/bin/time -l`;
with R it includes R's process).

Inputs:
- **GIAB HG008** (public), the pipeline's staged inputs and command,
  linked as they were staged.
- **A private whole-genome pair**, renamed copies (TUMOR / NORMAL) of the all-Java run's inputs (ids read from
  file names, never printed).
- **Synthetic inputs** (made from HG008's files): a cohort of 6 tumor/normal pairs with
  counts scaled at random, some lines dropped (so the samples' BQR contexts and GC buckets differ) and one
  ESVEE file missing; 37 edge directories (one file changed or left out each); threshold override files;
  cohort files (tumor rows only, header only, a null value).

**Against the pipeline's own outputs** (`qsee:1.1-beta.1`): the jar's three data files equal the pipeline's
for HG008 (byte for byte, gzip included) and for the private pair (after decompression; the pipeline's
copies were recompressed by the renaming). So 1.0 and 1.1-beta.1 write the same files on these inputs.

107 cases, all identical in Java mode (exit statuses, every data file byte for byte, gzip bytes included). Cases whose names start with `edge_` or `missing_` run HG008's files with one changed or left out; `missing_*_allow` with `-allow_missing_input`, `missing_*` without.

| Case | What | Java s, MB | Rust s, MB | Exit (J/R) | Result |
|---|---|---|---|---|---|
| hg008 | pipeline command, R report included | 12.26, 1349 | 9.02, 324 | 0/0 | identical |
| t0 | pipeline command, R report included | 14.15, 1438 | 9.63, 342 | 0/0 | identical |
| hg008_skip_plot | pipeline command, `-skip_plot` | 2.92, 1359 | 0.21, 34 | 0/0 | identical |
| t0_skip_plot | pipeline command, `-skip_plot` | 4.22, 1433 | 0.27, 38 | 0/0 | identical |
| hg008_tumor_only |  | 2.23, 1049 | 0.16, 31 | 0/0 | identical |
| t0_tumor_only |  | 3.77, 1391 | 0.19, 20 | 0/0 | identical |
| hg008_targeted |  | 2.90, 1337 | 0.22, 38 | 0/0 | identical |
| t0_targeted |  | 4.22, 1440 | 0.26, 33 | 0/0 | identical |
| hg008_sbx |  | 2.90, 1354 | 0.22, 39 | 0/0 | identical |
| hg008_ultima |  | 2.91, 1402 | 0.22, 34 | 0/0 | identical |
| hg008_bad_sequencing_type |  | 0.36, 94 | 0.01, 2 | 1/1 | identical |
| hg008_skip_categories |  | 2.84, 1299 | 0.20, 33 | 0/0 | identical |
| hg008_skip_summary |  | 0.36, 96 | 0.02, 2 | 1/1 | identical |
| hg008_skip_unknown |  | 0.34, 95 | 0.02, 2 | 1/1 | identical |
| hg008_lod_details | `-write_lod_details -output_id -log_debug -threads 2` | 3.03, 1408 | 0.31, 36 | 0/0 | identical |
| t0_lod_details |  | 4.35, 1433 | 0.36, 33 | 0/0 | identical |
| hg008_sample_data_dir | all inputs in one `-sample_data_dir` | 2.82, 1365 | 0.18, 24 | 0/0 | identical |
| hg008_no_panel | no driver gene panel | 0.58, 243 | 0.06, 14 | 0/0 | identical |
| hg008_no_cohort |  | 2.80, 1348 | 0.18, 24 | 0/0 | identical |
| hg008_bad_threads |  | 0.34, 93 | 0.02, 2 | 1/1 | identical |
| hg008_log_bad |  | 0.34, 90 | 0.02, 2 | 1/1 | identical |
| hg008_no_samples |  | 0.32, 89 | 0.01, 2 | 1/1 | identical |
| hg008_reference_only |  | 0.32, 90 | 0.01, 2 | 1/1 | identical |
| hg008_count_mismatch |  | 0.32, 91 | 0.01, 2 | 1/1 | identical |
| hg008_overrides |  | 2.89, 1358 | 0.21, 36 | 0/0 | identical |
| hg008_overrides_targeted |  | 2.89, 1365 | 0.22, 39 | 0/0 | identical |
| hg008_overrides_duplicate |  | 2.93, 1363 | 0.22, 36 | 0/0 | identical |
| hg008_overrides_unsupported |  | 0.35, 94 | 0.02, 2 | 1/1 | identical |
| hg008_overrides_purple |  | 2.89, 1351 | 0.21, 34 | 0/0 | identical |
| hg008_overrides_infinite |  | 0.61, 208 | 0.06, 28 | 1/1 | identical |
| hg008_overrides_java_test |  | 0.59, 218 | 0.06, 27 | 1/1 | identical |
| hg008_overrides_bad_operator |  | 0.35, 94 | 0.02, 2 | 1/1 | identical |
| hg008_cohort_tumor_rows |  | 0.58, 207 | 0.05, 19 | 1/1 | identical |
| hg008_cohort_header_only |  | 0.52, 189 | 0.03, 5 | 1/1 | identical |
| hg008_cohort_null |  | 0.52, 189 | 0.03, 5 | 1/1 | identical |
| hg008_tumor_only_cohort_tumor_rows |  | 2.20, 1031 | 0.14, 27 | 0/0 | identical |
| cohort_pairs | 6 synthetic pairs, `-sample_id_file`, wildcard dirs | 13.84, 1846 | 1.03, 77 | 0/0 | identical |
| cohort_pairs_t2 | the same at `-threads 2` | 13.76, 2081 | 0.98, 77 | 0/0 | identical |
| cohort_tumor_only |  | 10.97, 1802 | 1.21, 61 | 0/0 | identical |
| cohort_unsorted | sample file in reverse order, ReferenceId first | 14.02, 2089 | 1.03, 84 | 0/0 | identical |
| cohort_strings | 3 pairs by `-tumor`/`-reference` lists | 7.37, 1519 | 0.54, 60 | 0/0 | identical |
| cohort_missing_no_allow |  | 0.62, 211 | 0.04, 6 | 1/1 | identical |
| train_default | CohortPercentilesTrainer, 101 percentiles | 1.33, 340 | 1.04, 16 | 0/0 | identical |
| train_features | `-write_cohort_features -num_percentiles 5 -threads 2` | 0.90, 278 | 0.17, 10 | 0/0 | identical |
| train_increment |  | 0.64, 233 | 0.08, 5 | 0/0 | identical |
| train_increment_tenth |  | 0.37, 94 | 0.02, 2 | 1/1 | identical |
| train_bad_increment |  | 0.38, 96 | 0.02, 2 | 1/1 | identical |
| train_one_percentile |  | 0.74, 231 | 0.06, 6 | 0/0 | identical |
| train_missing_no_allow |  | 0.70, 237 | 0.05, 9 | 1/1 | identical |
| train_pair | trainer on the HG008 pair | 0.88, 283 | 0.77, 11 | 0/0 | identical |
| thresholds_default |  | 0.37, 92 | 0.01, 2 | 0/0 | identical |
| thresholds_targeted |  | 0.37, 92 | 0.01, 2 | 0/0 | identical |
| thresholds_all |  | 0.38, 94 | 0.01, 2 | 0/0 | identical |
| thresholds_all_targeted |  | 0.39, 92 | 0.01, 2 | 0/0 | identical |
| thresholds_cwd |  | 0.38, 92 | 0.01, 2 | 0/0 | identical |
| lod_entry | GeneLodCalculator entry point | 2.80, 1371 | 0.17, 21 | 0/0 | identical |
| lod_entry_details |  | 4.30, 1600 | 0.32, 22 | 0/0 | identical |
| edge_header_only_coverage |  | 0.41, 100 | 0.01, 2 | 1/1 | identical |
| edge_header_only_frag_length |  | 0.50, 169 | 0.02, 4 | 1/1 | identical |
| edge_header_only_dup_freq |  | 0.47, 160 | 0.01, 3 | 1/1 | identical |
| edge_header_only_bqr |  | 2.87, 1309 | 0.18, 24 | 0/0 | identical |
| edge_header_only_gene_coverage |  | 1.24, 533 | 0.08, 18 | 0/0 | identical |
| edge_empty_gc_median |  | 2.91, 1352 | 0.19, 28 | 0/0 | identical |
| edge_gc_median_zero |  | 2.92, 1362 | 0.18, 26 | 0/0 | identical |
| edge_disc_stats_zero_total |  | 2.90, 1368 | 0.18, 25 | 0/0 | identical |
| edge_disc_stats_short |  | 0.48, 160 | 0.01, 4 | 1/1 | identical |
| edge_disc_stats_new_name |  | 2.90, 1351 | 0.18, 27 | 0/0 | identical |
| edge_bqr_high_qual |  | 0.48, 156 | 0.02, 4 | 1/1 | identical |
| edge_bqr_lower_case |  | 2.84, 1357 | 0.18, 28 | 0/0 | identical |
| edge_bqr_no_read_type |  | 2.85, 1359 | 0.17, 28 | 0/0 | identical |
| edge_bqr_dual_single |  | 2.85, 1353 | 0.17, 24 | 0/0 | identical |
| edge_jitter_bad_unit |  | 0.50, 172 | 0.02, 4 | 1/1 | identical |
| edge_jitter_no_consensus |  | 2.86, 1363 | 0.18, 21 | 0/0 | identical |
| edge_jitter_old_name |  | 2.80, 1346 | 0.18, 23 | 0/0 | identical |
| edge_jitter_zero_counts |  | 2.82, 1354 | 0.18, 21 | 0/0 | identical |
| edge_male |  | 2.89, 1376 | 0.18, 26 | 0/0 | identical |
| edge_purity_low |  | 2.86, 1350 | 0.18, 23 | 0/0 | identical |
| edge_no_tumor |  | 2.89, 1369 | 0.17, 24 | 0/0 | identical |
| edge_warn_contamination |  | 2.86, 1349 | 0.18, 26 | 0/0 | identical |
| edge_fail_contamination |  | 2.90, 1359 | 0.17, 26 | 0/0 | identical |
| edge_tinc |  | 2.90, 1368 | 0.17, 23 | 0/0 | identical |
| edge_old_qc_names |  | 2.85, 1355 | 0.18, 24 | 0/0 | identical |
| edge_bad_purity |  | 0.41, 100 | 0.01, 2 | 1/1 | identical |
| edge_gene_on_chrM |  | 2.85, 1365 | 0.18, 24 | 0/0 | identical |
| edge_summary_zero_reads |  | 2.87, 1356 | 0.18, 25 | 0/0 | identical |
| edge_not_utf8 |  | 2.88, 1340 | 0.17, 25 | 0/0 | identical |
| missing_no_bam_summary_allow |  | 2.89, 1371 | 0.17, 22 | 0/0 | identical |
| missing_no_bam_summary |  | 0.37, 98 | 0.01, 2 | 1/1 | identical |
| missing_no_esvee_cobalt_allow |  | 2.89, 1363 | 0.18, 23 | 0/0 | identical |
| missing_no_esvee_cobalt |  | 0.48, 161 | 0.02, 4 | 1/1 | identical |
| missing_no_bqr_allow |  | 0.45, 161 | 0.02, 4 | 1/1 | identical |
| missing_no_bqr |  | 0.45, 162 | 0.02, 4 | 1/1 | identical |
| missing_no_exon_coverage_allow |  | 1.23, 605 | 0.08, 16 | 0/0 | identical |
| missing_no_exon_coverage |  | 0.54, 205 | 0.03, 5 | 1/1 | identical |
| missing_no_purple_allow |  | 1.22, 534 | 0.08, 18 | 0/0 | identical |
| missing_no_purple |  | 0.37, 99 | 0.01, 2 | 1/1 | identical |
| missing_no_purple_qc_allow |  | 2.80, 1355 | 0.18, 25 | 0/0 | identical |
| missing_no_purple_qc |  | 0.38, 97 | 0.01, 2 | 1/1 | identical |
| missing_no_redux_normal_allow |  | 0.49, 196 | 0.02, 4 | 1/1 | identical |
| missing_no_redux_normal |  | 0.51, 198 | 0.02, 4 | 1/1 | identical |
| missing_only_purple_and_summary_allow |  | 0.39, 100 | 0.01, 2 | 1/1 | identical |
| missing_only_purple_and_summary |  | 0.37, 99 | 0.01, 2 | 1/1 | identical |
| missing_not_utf8_allow |  | 2.87, 1338 | 0.17, 21 | 0/0 | identical |
| missing_not_utf8 |  | 0.40, 101 | 0.01, 2 | 1/1 | identical |
| no_dirs_allow |  | 0.41, 120 | 0.01, 3 | 1/1 | identical |
| no_dirs | inputs in the working directory, no directory options | 0.39, 98 | 0.01, 2 | 1/1 | identical |
| wildcard_dollar | `-tumor S$1` with a wildcard directory | 0.34, 91 | 0.01, 2 | 1/1 | identical |

### The corrected (default) mode against Java

The regression in its default mode, summarised by a script. Apart from QS4's one display change (the
CONTAMINATION threshold shown as `>=`, on every run with PURPLE features), the default mode writes Java's
files except in these cases, each a finding:

| Case | Default mode | Finding |
|---|---|---|
| missing_no_bqr_allow, missing_no_redux_normal_allow, missing_only_purple_and_summary_allow | exit 0 and the three files (BQR features left out) where Java exits 1 | QS1 |
| thresholds_default, _targeted, _all, _all_targeted, _cwd | 10 threshold lines (the other mode's) | QS2 |
| hg008_cohort_tumor_rows, hg008_cohort_header_only | exit 0, the three files, warnings | QS3 |
| every case with PURPLE features | the CONTAMINATION condition `>=0.1` / `>=10%` | QS4 |
| edge_warn_contamination | CONTAMINATION `WARN >=0.02` where Java has `PASS >0.1` | QS5 |
| hg008_overrides_infinite, hg008_overrides_java_test | exit 0, the three files (Java: exit 1, a vis data file of 2 lines) | QS6 |
| no_dirs, no_dirs_allow | exit 0, the inputs found in the working directory | QS7 (and QS1) |
| hg008_overrides_purple | PURITY and CONTAMINATION keep PURPLE's conditions | QS9 |
| wildcard_dollar | exit 0 (inputs missing, allowed) where Java fails on the replacement | QS11 |
| train_increment_tenth | the percentiles file (1,001 percentiles) where Java exits 1 | QS12 |

Every other case gives Java's files in the default mode too, QS4's line apart.

## The report, drawn natively

`<tumor>.qsee.vis.report.pdf` and (single patient) `.png` are drawn by `src/plot/` with hmf-plot-rs (`ggp`,
`gtable`, `patchwork`, `rdev`, `pdf_figure`), as the repository's `plot_qc.R` draws them with R 4.5.3, ggplot2
4.0.3 and patchwork 1.3.2 in the qsee 1.0 image (whose own jar has an older `plot_qc.R`; the references run the
repository's script there, as the pipeline's qsee 1.1-beta.1 does). Multi-sample mode writes one PDF per tumor
(no PNG), and `-merge_plots` one `multisample.qsee.vis.report.pdf` with the per-sample PDFs deleted, as Java.

How it was checked (2026-11-27, the port in Java mode): R drew each case in the image (docker, amd64 emulated),
and again with every grid viewport's device box dumped, for the PNG (300 dpi) and the PDF (points); the port wrote
its grobs' boxes and its label boxes. Layout: boxes paired by name (gtable and grid's own levels dropped), then
every box without a name partner looked for by position among all of the other side's. PNG: pixels compared.
PDF: pages and page size, `pdftotext -layout`, every word's box against R's, the page at 100 dpi.

| case | R in the image: s, MB | port: s, MB | PNG pixels differing by >10% / >25% / >50% (of 15.3 M) | layout elements, largest difference (PNG and PDF) | PDF: text, words at R's position, pixels >25% at 100 dpi |
|---|---|---|---|---|---|
| hg008 | 18.58, 298 | 0.31, 186 | 59,617 / 24,900 / 3,073 | 728, 0.0000 px | identical, 569 (0.000 pt), 2 |
| t0 (private pair) | 18.88, 299 | 0.31, 183 | 60,368 / 26,662 / 2,938 | 737, 0.0000 | identical, 609, 13 |
| hg008_tumor_only | 17.89, 320 | 0.28, 179 | 64,014 / 26,689 / 3,151 | 728, 0.0000 | identical, 556, 3 |
| t0_tumor_only | 14.53, 336 | 0.24, 174 | 61,947 / 26,388 / 2,973 | 737, 0.0000 | identical, 602, 15 |
| hg008_no_cohort | 16.01, 303 | 0.25, 175 | 59,796 / 24,986 / 3,088 | 728, 0.0000 | identical, 570, 2 |
| edge_purity_low | 14.84, 304 | 0.26, 176 | 58,072 / 25,059 / 3,047 | 737, 0.0000 | identical, 594, 4 |
| edge_fail_contamination | 15.35, 303 | 0.25, 175 | 60,070 / 25,082 / 2,993 | 728, 0.0000 | identical, 575, 3 |
| edge_disc_stats_zero_total | 14.75, 336 | 0.25, 176 | 60,259 / 25,889 / 3,253 | 728, 0.0000 | identical, 549, 2 |
| missing_no_purple_allow (QS14, Java mode) | exit 1, 4.44, 211 | exit 1, 0.02, 6 | no PNG (as R) | | blank page as R's, 0 |
| missing_no_purple_fixed (QS14 corrected, R patched) | 13.27, 310 | 0.24, 175 | 50,802 / 21,302 / 2,830 | 728, 0.0000 | identical, 517, 1 |
| missing_no_esvee_cobalt_allow | 14.67, 336 | 0.26, 180 | 60,148 / 25,832 / 3,244 | 728, 0.0000 | identical, 549, 2 |
| missing_no_bam_summary_allow | 14.83, 303 | 0.25, 175 | 59,581 / 24,946 / 3,088 | 728, 0.0000 | identical, 564, 2 |
| edge_header_only_bqr | 14.78, 332 | 0.24, 175 | 59,713 / 24,977 / 3,088 | 728, 0.0000 | identical, 570, 2 |
| hg008_targeted | 19.07, 298 | 0.31, 180 | 63,137 / 27,109 / 3,168 | 728, 0.0000 | identical, 584, 9 |

- Every port box has an R viewport at the same place (0.01 px). The R viewports without a port box draw nothing:
  the legend patch's empty title, subtitle, caption and tag cells, legend.box.background, the QC-status keys'
  zero-size anchors, and (three cases) the second facet's axis viewports, which share three edges with the cell
  the port records.
- Multi-sample mode (6 tumor/normal pairs, private): 0.37 s and 44 MB for the 6 PDFs (R: about 11 s and 225 MB
  per sample); samples 1 and 6 against R: text identical, every word at R's position, 5 and 2 pixels differ at
  100 dpi; `-merge_plots`: 6 pages, page i sample i's report (0.46 s, 35 MB).
- Visible differences (PNG; none at normal size): text is rasterised by hmf-plot-rs with whole-pixel glyph
  positions (Cairo positions glyphs to a fraction of a pixel): horizontal text lands at R's position, the bold
  title 1 px lower (BUI's synthetic bold), and a rotated axis title drifts by up to 2 px along its length; panels
  clip to whole pixels, so a panel border differs in its outermost pixel line (13% of the pixels differing by
  >25% on hg008, none by >50%); anti-aliasing is BUI's. QS15's cut strip label is as in R.
- Fixed on 2026-11-27 after the first comparison: the summary table's extra nesting level that `free()` adds
  (`merge()`), the geom_label boxes (half a descent too low), the QC-status keys' border width, the legend key
  backgrounds, the PDF's dash lengths for lines thinner than lwd 1, rectangles in the PDF written as R writes them,
  and empty (zero-height) panels, which drew a line (hmf-plot-rs README).
- Determinism: a second run gives the same bytes; built on Linux arm64 (`rust:1.97-bookworm`) the HG008 PNG and
  PDF have the same bytes as on macOS arm64.
- Speed and memory against native R on the same Mac (one thread, three runs): 0.31 s and 181 MB for HG008's PDF
  and PNG against 8.8-9.1 s and 321 MB.
- ORANGE: with the port's report in `qsee/`, ORANGE's JSON has `qSeePlot` as Java's (byte-identical to Java's
  JSON on HG008) and its PDF has the qsee page (9 pages as Java's); without it the JSON had null there and the
  PDF 8 pages.

## Findings

Each corrected by default behind `compat::java("QSn")`, `-java_compat` gives Java's behaviour; QS8 and QS10
have no switch (see each).

- **QS1 (BaseQualRecalibrationPrep / BqrFile.read): a missing BQR file ends the run, even with
  `-allow_missing_input`.**
  - **What Java does:** BqrFile.read catches every exception, logs it and returns null; loadSnvBqrRecords
    then streams the null list, and the NullPointerException is not the IOException CategoryPrepTask
    treats as a missing input. Its catch-all logs "Failed to run prep" and calls System.exit(1).
  - **Effect:** a sample without its REDUX BQR file cannot be run at all; every other input is optional with
    `-allow_missing_input`. In tumor/normal mode a missing normal BQR file ends the tumor's run too.
  - **Evidence:** `missing_no_bqr_allow`, `missing_no_redux_normal_allow`, `missing_only_purple_and_summary_allow`
    and `no_dirs_allow`: Java exit 1 and no files; the default mode writes the three files with the BQR
    features left out and a warning.
  - **Corrected:** a BQR file that cannot be read is a missing input; one that cannot be parsed still ends
    the run (as Java).
- **QS2 (DefaultThresholdsWriter.main): `-targeted_mode` is inverted.**
  - **What Java does:** `boolean targetedMode = !configBuilder.hasFlag(TARGETED_MODE)`.
  - **Effect:** without the flag the file has the targeted thresholds, with it the WGS ones (the README says
    to add the flag for the targeted thresholds). 10 lines differ (duplicate reads, dual-strand reads and
    the tumor coverage thresholds).
  - **Evidence:** the five `thresholds_*` cases.
  - **Corrected:** the flag selects the targeted thresholds. `-no_overridable_only` (inverted in name and
    meaning, as its description says) is left as it is.
  - **Fixed upstream** in bdceb2302a (`targetedMode = configBuilder.hasFlag(TARGETED_MODE)`): since the port follows
    883ab11bba both modes take the flag as the correction did, and the switch is gone.
- **QS3 (CohortPercentiles.getFeaturePercentiles): a cohort file without rows of a sample type crashes the
  run.**
  - **What Java does:** `mCohortData.get(sampleType).get(featureKey)`; with no row of that type the first
    `get` is null, and the NullPointerException escapes QseeApplication.run.
  - **Effect:** a cohort percentiles file trained on tumors only (or a header-only file) cannot be used for
    a tumor/normal run (`-cohort_percentiles_file` is meant to be optional and only adds warnings and the
    plot's ranges). Java exits 1 before writing anything.
  - **Evidence:** `hg008_cohort_tumor_rows`, `hg008_cohort_header_only` (exit 1); the tumor-only run with
    the same file works (`hg008_tumor_only_cohort_tumor_rows`).
  - **Corrected:** a missing sample type is a missing percentile: a warning per feature.
- **QS4 (ThresholdRegistry.setDefaults): the contamination thresholds are shown as `>` where PURPLE uses
  `>=`.**
  - **What Java does:** CONTAMINATION's WARN and FAIL thresholds are registered with GREATER_THAN, but
    PurpleQCStatus sets FAIL_CONTAMINATION at `contamination >= 0.1` and WARN_CONTAMINATION at `>= 0.02`
    (the TINC thresholds follow PURPLE: `>` 0.02 and `>=` 0.25).
  - **Effect:** only the text: FailCondition `>0.1` and the vis data's `QcThreshold:>10%`, on every run with
    PURPLE features. A sample at exactly 10% shows FAIL with the condition `>0.1` (`edge_fail_contamination`).
  - **Corrected:** `>=`.
- **QS5 (SummaryTablePurplePrep.createFeatures): PURPLE's WARN_CONTAMINATION is shown as PASS.**
  - **What Java does:** CONTAMINATION takes the mapping of FAIL_CONTAMINATION only. TINC, the other
    two-level PURPLE check, takes FAIL, else WARN, else PASS (getTincQcStatus).
  - **Effect:** a sample PURPLE warns for (contamination from 2% to 10%) passes in qsee's status file and
    report, with the condition `>0.1`.
  - **Evidence:** `edge_warn_contamination` (contamination 0.02, PURPLE status WARN_CONTAMINATION): Java
    `PASS >0.1`, default mode `WARN >=0.02`.
  - **Corrected:** as TINC.
- **QS6 (QcStatus.formDiplayString): an infinite threshold ends the run while the files are written.**
  - **What Java does:** BigDecimal.valueOf(-Infinity) throws NumberFormatException. ThresholdOverridesFile
    reads `-Infinity` (or `Infinity`) without complaint, and the exception comes when the vis data file
    shows the threshold.
  - **Effect:** exit 1, a vis data file of the lines before the feature (2 of 3,847 for HG008: the header and
    MEAN_COVERAGE), no status or LoD file. Java's own test resource (`threshold_overrides.tsv`, MAPPED_PROPORTION FAIL < -Infinity,
    the way to switch a check off) does it.
  - **Evidence:** `hg008_overrides_infinite`, `hg008_overrides_java_test`.
  - **Corrected:** shown as `<-Inf` (`<-Inf%` for percentages).
- **QS7 (PurplePurity, PurpleQCFile, CobaltGcMedianFile, EsveeDiscordantStats, BqrFile,
  JitterCountsTableFile generateFilename): without a tool directory the files are looked for at the file
  system root.**
  - **What Java does:** these join `basePath + File.separator + file`; with neither `-<tool>_dir` nor
    `-sample_data_dir` the base is "" and the path `/<sample>.purple.purity.tsv`. The BamMetrics and duplicate
    frequency names use checkAddDirSeparator, which leaves "" alone, so those are read from the working
    directory.
  - **Effect:** with every input in the working directory and no directory options, Java finds the BamTools
    files but none of PURPLE's, COBALT's, ESVEE's or REDUX's BQR and jitter files (exit 1 without
    `-allow_missing_input`, and with it exit 1 through QS1).
  - **Evidence:** `no_dirs`, `no_dirs_allow` (run in a directory holding the inputs).
  - **Corrected:** an empty directory is the working directory for every file.
- **QS8 (FeaturePrep.prepCohortCategory / FeatureMatrix.addRow): with `-threads` above 1, multi-sample
  output order depends on thread timing.**
  - **What Java does:** the samples of a category run on a thread pool; each adds its features to the
    matrix when it finishes, and a key's column is placed when first added. sortFeatureKeys then sorts, but
    FeatureKey.compareTo compares BQR features by their quality bin only, so the BQR_PER_SNV96_CONTEXT keys
    (all one bin) keep the order in which the samples happened to finish.
  - **Effect:** where the samples' BQR contexts differ, the vis data's row order changes from run to run (the
    trainer's percentiles and features files likewise). The values are the same.
  - **Evidence:** six Java runs of `cohort_pairs_t2`: three wrote the order of the
    1-thread run, three another order (18 BQR_PER_SNV96_CONTEXT lines moved, same lines).
  - **The port** adds each category's samples in sample order whatever `-threads` is, which is Java's
    1-thread order: the defined Java reference, in both modes (no switch).
- **QS9 (ThresholdOverridesFile.read): overrides of thresholds PURPLE decides change only the threshold
  shown.**
  - **What Java does:** the PURPLE features take their status from PURPLE's QC statuses and only the shown
    threshold from the registry; an override of, say, TUMOR PURITY WARN is accepted (DefaultThresholdsWriter
    marks these "Overridable false" but the reader does not check).
  - **Effect:** a purity of 0.25 with an override `< 0.3` passes, shown with the condition `<0.3`.
  - **Evidence:** `hg008_overrides_purple`: Java shows `<0.3` and `>0.05` for PURITY and CONTAMINATION while
    the statuses are PURPLE's.
  - **Corrected:** such an override is ignored with a warning, so the condition shown is PURPLE's.
- **QS10 (SummaryTablePurplePrep): FAIL_NO_TUMOR and WARN_GENDER_MISMATCH never reach the QC status.**
  - **What Java does:** both PURPLE statuses have thresholds registered (null operator, no value), but no
    feature takes them. A sample PURPLE calls NO_TUMOR shows PURITY PASS unless WARN_LOW_PURITY is also
    set, which PURPLE does not set for NO_TUMOR.
  - **Effect:** the report's overall status can be PASS for a sample PURPLE failed (`edge_no_tumor`).
  - **Not corrected:** which feature should carry them (purity? a new row?) is not clear from the code, so
    Java's behaviour is kept; listed for the authors.
- **QS11 (ConfigUtils.convertWildcardSamplePath): a `$` or `\` in a sample id breaks a wildcard directory.**
  - **What Java does:** `samplePath.replaceAll("\\*", sampleId)` takes the id as a regex replacement: `$1`
    is a group reference (IndexOutOfBoundsException "No group 1"), `$0` the `*` itself, `\x` an escape.
  - **Evidence:** `wildcard_dollar` (`-tumor S$1` with `-bam_metrics_tumor_dir .../*/`): Java exit 1.
  - **Corrected:** the id is inserted literally.
- **QS12 (PercentileTransformer.withIncrement): `-percentile_increment` rejects 0.1.**
  - **What Java does:** `100 % increment != 0` in floating point; 100 % 0.1 is 0.09999999999999445, so 0.1,
    0.2 and other increments that are not exact in binary are rejected, though the error message gives 0.1
    as an example of a valid one.
  - **Evidence:** `train_increment_tenth`: Java exit 1 "Percentile increment should evenly divide 100".
  - **Corrected:** the increment must divide 100 to within 1e-7; 0.1 gives 1,001 percentiles (as
    `-num_percentiles 1001`).
- **QS13 (plot_qc.R, coverage and fragment length panels): the mirrored y axis loses its decimals.**
  - **What R does:** the normal sample is drawn below zero and the axis labelled
    `scales::label_percent()(abs(x))`. The breaks carry rounding noise (`abs(-0.025)` and `abs(0.025)` differ
    in the last bits), so `scales::precision` sees a gap below sqrt(machine eps) between them and labels at
    accuracy 1: 2.5% on both sides becomes "2%" above zero and "3%" below, and an axis of 0.25% steps reads
    "0%" at every break.
  - **Evidence:** HG008's report (R's in the qsee image, as the pipeline's): the coverage axis, in steps of
    2.5%, reads 2%, 0%, 3%, 5%, and the fragment length axis 0%, 0%, 0%; the port in Java mode draws the same
    labels.
  - **Corrected:** the absolute values are taken at 12 significant digits before labelling, so the labels
    have the accuracy the breaks need: on HG008 2.5%, 0.0%, 2.5%, 5.0% and 0.25%, 0.00%, 0.25%. The breaks
    are unchanged; the wider labels narrow those two panels by their extra width.
- **QS14 (plot_qc.R with ggplot2 4.0.3): without PURPLE's outputs the report is not drawn and qsee fails.**
  - **What R and Java do:** with `-allow_missing_input` and no PURPLE files, the summary table's two LOG10
    sub-tables (Copy number, Mutational burden) have no rows. Their dummy limit points inherit the empty data,
    so the log scale is never trained; ggplot2 4.0.3's `guide_axis_logticks` then gets no parameters, and
    under `coord_flip` its transform calls `setNames(NULL, ...)`: "Error in names(object) <- nm : attempt to
    set an attribute on NULL". R leaves a blank one-page PDF and no PNG; QseePlotTask throws a
    RuntimeException, so qsee exits 1 after writing its data files.
  - **Evidence:** case missing_no_purple_allow (the Java baseline's vis data, R in the qsee image): R exit 1;
    the call stack ends in `GuideAxisLogticks$transform` -> `CoordFlip$transform` -> `flip_axis_labels` ->
    `setNames`. The port in Java mode: exit 1, the same blank page (no pixel differs at 100 dpi), no PNG.
  - **Corrected:** the report is drawn as R draws it once the guide skips an empty scale (as `guide_axis`
    does; `scripts/ggplot2_logticks_fix.R` outside this repository applies that to R for the reference): the
    empty sub-tables keep their titles and zero-height panels (R's y axis there, trained only by the div lines
    at `seq(1.5, 0)`, prints 0.50 to 1.50 on top of each other, and the port draws them so). Case
    missing_no_purple_fixed: the layout equals patched R's, pdftotext identical, every word at R's position.

Kept as R draws it, recorded:
- **QS15 (plot_qc.R): the MS indel panels' last strip label is cut off.** The four repeat-type strips of
  panels J and K share the panel width, and ">=3bp repeat" is wider than its strip at 8.4 pt, so R clips both
  of its ends in the PNG (DejaVu Sans; Helvetica in the PDF is narrower). Fixing it needs a layout choice
  (a smaller strip font, shorter labels or wider panels) that belongs to the script; the port clips as R.
- **QS16 (ggplot2 4.0.3's draw_key_label, in the report's QC-status legend): a size in mm used as points.** The
  key is raised by `font_descent(size = data$size)`, the label size in mm passed where a font size in points is
  expected (the key's text is set at `size * .pt`), so the raise is the descent of a font about a third of the
  size: 3.7 px at 300 dpi in R's viewports, where the label's own descent is about 9 px. The intent is unclear:
  with 4.0.3's labelGrob, which centres the box on its point, no raise at all would centre the box in its key
  cell, so the slip brings it closer to the centre than consistent units would. Cosmetic; the port draws the keys
  where R does (hmf-plot-rs `ggp::geom::key_label`, tested against R's box).

- **QS17 (plot_qc.R at 883ab11bba with ORANGE at the same commit): the report's PNG is renamed, and ORANGE loses
  its qsee page.**
  - **What R does:** the report has two pages; for a single patient each is written as
    `<tumor>.qsee.vis.report.page_<n>.png` (`sub("pdf$", sprintf("page_%s.png", i), PLOT_PATH)`), and
    `<tumor>.qsee.vis.report.png` is no longer written.
  - **Effect:** ORANGE (and orange-rs) embeds `<qsee_dir>/<tumor>.qsee.vis.report.png`
    (hmf-common MiscToolFiles.QSEE_PLOT_FILE_ID, unchanged at 883ab11bba); OrangeAlgo leaves qSeePlot null when the
    file is missing, so the JSON's `plots.qSeePlot` is null and the PDF's Quality control chapter says "Plot not
    available". The commit does not update ORANGE, and oncoanalyser's QSEE module publishes the whole `qsee/`
    folder, so nothing else picks the new names up.
  - **Evidence:** the R references of this review (`work/up_pave_qsee/refs/qsee/<case>/`: page_1.png and page_2.png
    only); OrangeAlgo.java lines 440-448 and QualityControlChapter.java at 883ab11bba.
  - **Corrected:** page 1 (the summary page, what the single page was) is also written as
    `<tumor>.qsee.vis.report.png` in the default mode; `-java_compat` writes the two page files only.

Remarks, not findings:
- The report's legend is built while R's `pdf(NULL)` device is open, so in the PNG too its keys and texts
  are sized with Helvetica's metrics and drawn in DejaVu Sans (slightly wider); the port does the same.
- `plot_qc.R` in the qsee 1.0 biocontainer's jar is older than this repository's (and than the pipeline's
  qsee 1.1-beta.1); the port follows the repository's.
- The discordant fragment features divide fragment counts by total reads (ESVEE notes its counts are in
  fragment terms, "hence see doubling in rate calc"), so they are half the share of reads; the cohort
  percentiles use the same definition, so the comparison is consistent and nothing is changed.
- `-output_id` is registered with the description "Output directory" (QseeConfig), as `-help` shows.
- `-num_percentiles 1` gives one percentile NaN (100 / 0 steps); the file has a column `Pct_NaN`.
- `compare_cohorts.R` is in the jar but no Java class runs it.

## Following upstream (883ab11bba, reviewed 2026-11-28)

The port started from ca04e6c39d. Upstream master 883ab11bba changes qsee twice: bdceb2302a (exon coverage QC) and
883ab11bba (read strands in the BQR features and report). hmf-common's changes in that range (PipelineToolDirectories,
VirusInterpreterDataLoader, from f5f660cc9f) are code qsee does not use. The decisions, with reasons and evidence,
are in the review's NOTES (work/up_pave_qsee/port/NOTES.md) and UPSTREAM.md; in short:

- **bdceb2302a, accepted:** the gene LoD file becomes `<tumor>.qsee.exon_coverage_qc.tsv.gz` (LodVaf, the sample VAF
  detected 95% of the time, floored at 0.01; PASS / WARN against the new EXON_VARIANT_LOD threshold, tumor 0.1,
  targeted 0.05, normal 0.3; EXCLUDED); "Low coverage exons" (the share of WARN exons, WARN above 10%) joins the
  summary table; the missed-variant likelihood comes from the exon coverage file with excluded exons left out
  (unchanged values without exclusions: BamTools' gene counts are its exon counts summed); `-lod_min_alt_reads`,
  `-exon_exclusions_file` and the ExonExclusionsFileBuilder entry point; DefaultThresholdsWriter's `-targeted_mode`
  no longer inverted, which is QS2's correction (the switch is gone).
- **883ab11bba, accepted with a correction (QS17):** BQR features by read strand (the C/T-centred context, REVERSE
  when the read carried the purine; REDUX's BQR keys are read-oriented), C>T split into CpG and other per original
  quality, groups of fewer than 5 counts NaN (kept in the vis data, hidden in the report); the report on two pages
  with the strands as triangles. R writes the pages' PNGs as `.page_<n>.png`, which leaves ORANGE without its qsee
  page (QS17): the default mode also writes page 1 as `.png`.

**Inputs added** for the new code (`baseline/qsee/inputs/up/`, made by `work/up_pave_qsee/scripts/qsee/make_up_inputs.py`
from the public HG008 inputs and public gene names): exon exclusion files (TP53's exons, two of KRAS, a gene outside
the panel; header only; a rank that is not a number), threshold overrides of the new thresholds, an HG008 directory
whose tumor has 60 exons at 1500x and 20 at 0x, and one whose BQR file has counts set to 0 and 1 (NaN groups). The
synthetic cohort (S1-S6) already has exon coverage files.

**Cases and results.** The 107 cases of `work/qsee/port/scripts/regress.sh` and 20 new ones (`up_*`: exclusions on
HG008, the private pair and the cohort, `-lod_min_alt_reads` 5, 1 and not a number, overrides of the new thresholds
in WGS and targeted mode, the high-depth exons, the low BQR counts, ExonExclusionsFileBuilder on the cohort, with
options, with one tumor and without `-allow_missing_input`, ExonCoverageQcFileBuilder with exclusions, the trainer on
low BQR counts) were run with the 883ab11bba jar (`baseline/qsee/java_883ab11bba/<case>`; the two pipeline cases with
the new R script on the Mac's R) and with the port in Java mode (`baseline/qsee/rust_up/<case>`); script
`work/up_pave_qsee/scripts/qsee/regress.sh`. **126 of 127 identical** (exit statuses, every data file byte for byte,
gzip bytes included). The other, `cohort_pairs_t2` (6 pairs at `-threads 2`), is QS8: Java's 2-thread file has the
same lines in another order (now including NaN BQR lines of contexts some samples lack); the port's files are
byte-identical to Java's at `-threads 1`, the defined reference (`java_883ab11bba/cohort_pairs_t2_t1`,
`scripts/qsee/qs8_check.sh`). The default mode differs from Java exactly where the findings say (QS1, QS3-QS7, QS9,
QS11, QS12; QS2 no longer: `thresholds_*` are as Java), as before (`scripts/qsee/default_diff_summary.py`).

**Effect in the default mode** (the merged port, `rust_before_default`, against this one, `rust_up_default`;
`scripts/qsee/effect.py`, log `work/up_pave_qsee/logs/qsee_effect.txt`), on every case that runs:
- `gene_lod` is replaced by `exon_coverage_qc` (same exons, the same LoD where not floored; none of our WGS exons is
  floored); the statuses: HG008 tumor 3,478 PASS, 2,650 PASS -> WARN, 8 FAIL -> WARN, 3 WARN_HET_UNRELIABLE -> WARN;
  HG008 normal 5 PASS -> WARN, 1 WARN_HET_UNRELIABLE -> WARN; private pair tumor 1,647 PASS, 3,601 PASS -> WARN, 171
  FAIL -> WARN, 720 WARN_HET_UNRELIABLE -> WARN.
- the vis data and status files gain LOW_COVERAGE_EXONS per sample: HG008 tumor 43.3% (WARN), normal 0.6%; private
  pair tumor 73.2% (WARN), normal 0.1%; HG008 targeted tumor 86.8% (WARN, threshold 0.05). The report's title then
  lists "Low coverage exons (>10%)" for both tumors.
- BQR_PER_ORIG_QUAL 48 -> 112 lines and BQR_PER_SNV96_CONTEXT 384 -> 768 per pair (the new names); the pipeline's
  cohort percentiles file has none of them, so 880 more "Missing cohort percentiles" warnings per pair (915 against
  33) and no BQR cohort boxes in the report until a cohort file is trained on the new names.
- MISSED_VARIANT_LIKELIHOOD unchanged everywhere (values and order).
- `missing_no_exon_coverage` (no `-allow_missing_input`): the summary table now needs the exon file, so the run ends
  before writing the vis data and status files (as Java); the trainer's files gain the new features; the default
  thresholds file gains the 4 new thresholds.
- the report: two pages (and page 1 also as `.png`, QS17).

**The report.** R draws each case in the qsee image (`work/up_pave_qsee/scripts/plots/r_ref.sh`, the PNG per page and
the PDF; `r_layout.sh`, each page on a PDF and a 300 dpi PNG device with every viewport's box dumped), the port in
Java mode with its boxes (`compare.sh`, `plot_cases.sh`); side by side (R left) in `work/up_pave_qsee/compare/qsee_<case>/`
(private cases under `baseline/up_pave_qsee/compare/`).

| case | R in the image: exit, s, MB | port: exit, s, MB | PNG px >10/25/50% page 1; page 2 | layout elements, >0.01 px, largest: page 1; page 2 | PDF pages R/port, text, words (largest pt), px >25% at 100 dpi p1/p2 |
|---|---|---|---|---|---|
| hg008 | 0, 73.61, 354 | 0, 0.52, 205 | 52098 / 20532 / 2723; 44998 / 10064 / 1203 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 553 (0.000), 11/58 |
| t0 | 0, 79.70, 357 | 0, 0.51, 205 | 59268 / 25489 / 2216; 50011 / 13286 / 1094 | 697, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 616 (0.000), 18/60 |
| hg008_tumor_only | 0, 48.36, 350 | 0, 0.46, 201 | 50080 / 21148 / 2795; 42852 / 9702 / 1203 | 688, 0, 0.0000; 90, 0, 0.0000 | 2/2, identical, 525 (0.000), 8/33 |
| t0_tumor_only | 0, 22.61, 338 | 0, 0.43, 193 | 59201 / 24179 / 2131; 45520 / 11621 / 1101 | 697, 0, 0.0000; 90, 0, 0.0000 | 2/2, identical, 606 (0.000), 16/37 |
| hg008_no_cohort | 0, 22.49, 328 | 0, 0.46, 195 | 51979 / 20496 / 2730; 44998 / 10064 / 1203 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 555 (0.000), 11/58 |
| hg008_newcohort | 0, 40.38, 364 | 0, 0.54, 205 | 55184 / 24997 / 2795; 45019 / 10089 / 1203 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 559 (0.000), 15/41 |
| bqr_low_counts | 0, 42.57, 365 | 0, 0.53, 207 | 55317 / 25009 / 2795; 44925 / 10068 / 1203 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 555 (0.000), 12/40 |
| edge_header_only_bqr | 0, 20.64, 334 | 0, 0.43, 194 | 51741 / 20415 / 2730; 42744 / 9667 / 1203 | 688, 0, 0.0000; 90, 0, 0.0000 | 2/2, identical, 547 (0.000), 11/31 |
| edge_bqr_dual_single | 0, 22.62, 300 | 0, 0.45, 196 | 52424 / 21053 / 2733; 39151 / 8399 / 1203 | 713, 0, 0.0000; 132, 0, 0.0000 | 2/2, identical, 577 (0.000), 11/68 |
| hg008_targeted | 0, 69.28, 354 | 0, 0.51, 205 | 55572 / 22636 / 2807; 49411 / 12997 / 1305 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 583 (0.000), 20/58 |
| edge_fail_contamination | 0, 21.22, 328 | 0, 0.45, 194 | 52217 / 20530 / 2634; 45549 / 10399 / 1084 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 565 (0.000), 14/58 |
| missing_no_esvee_cobalt_allow | 0, 20.95, 330 | 0, 0.45, 194 | 53051 / 19192 / 1905; 44998 / 10064 / 1203 | 688, 0, 0.0000; 107, 0, 0.0000 | 2/2, identical, 540 (0.000), 9/58 |

The layout of every page of every case equals R's to 0.01 px (every element both name; the boxes left unpaired are
of the kinds listed in "The report, drawn natively": the legend patch's empty title, subtitle and caption cells,
legend.box.background, zero-size anchors); the PDFs have R's pages and size, `pdftotext -layout` identical and every
word at R's position. Pixel differences are those of the first round (glyph rasterisation, whole-pixel clips,
anti-aliasing); the triangles are drawn as R's (pch 24 and 25: the circle's area, apex at 1.5551 times its radius).
Checked by eye on HG008 and the low-count case: the same points, triangles, hidden NaN values and boxes.

## Speed and memory

The port is 14 (HG008) and 16 (private pair) times faster than Java on the pipeline's data files, with a
fortieth of the memory: 0.21 s and 34 MB against 2.92 s and 1,359 MB (HG008), 0.27 s and 38 MB against 4.22 s
and 1,433 MB (private pair), both with `-skip_plot`. With the R report (the pipeline's command) both take
about 9 s more, R's time: 9.0 s against 12.3 s (HG008), 9.6 s against 14.2 s; R's peak (about 330 MB) is then
the port's. In multi-sample mode (12 synthetic samples) 1.0 s against 13.8 s. (Those rows ran R. Since
2026-11-27 the port draws the report itself in about 0.1 s: 0.31 s and 181 MB for HG008 with the report, the
memory mostly the 5100 x 3000 canvas.)

Where the port's 0.19 s of CPU goes (HG008, from the log's time stamps):
- the tumor's ten categories: 7 ms; the normal's and the cohort file (1.1 MB gzip, 4,000 rows of 101
  percentiles): 32 ms;
- the gene limit of detection: about 0.1 s. Each of 7,230 exons of reportable genes bisects the missed
  variant likelihood 20 times, each a binomial CDF (a regularized beta) per depth range;
- writing the three files (gzip level 6): about 20 ms.

Java's 2.9 s includes the JVM's start and JIT compilation (4.6 s of CPU over its threads). Java uses
threads only for multi-sample mode and the trainer (`-threads`); the port too, adding the samples in order.
The trainer's time (1.0 s against Java's 1.3 s for 12 samples and 101 percentiles) goes into writing the
400,000 values of the percentiles file (hmf-common-rs's DecimalFormat port, then zlib); the trainer is not a
pipeline step.

## Shared code added to hmf-common-rs

Additions only; no existing behaviour changed. `cargo test --release` in hmf-common-rs: 122 pass (1 ignored), the 8 new ones included.

- `delim_file` (new): DelimFileReader (header as column names, a repeated name an error, `null`/`NULL`
  values an error unless asked for, Integer/Long/Double parsing), `JavaWriter` (FileWriterUtils.
  createBufferedWriter: plain, or GZIPOutputStream's bytes for `.gz`), `DelimRow` / `write_delim_file`
  (DelimFileWriter: unset values as "null", a repeated column an error, the file closed with the rows before
  a failing one). Tests `reader_like_java`, `writer_gzip_and_plain`.
- `redux_files` (new): ConsensusType, BqrFile.read (BqrRecord; case-insensitive columns, (byte) qualities),
  DuplicateFrequency.read, JitterCountsTableFile.read (JitterCountsTable, JitterTableRow, rows kept in order of
  units), with the file names. Tests `bqr_lines`, `duplicate_frequency_lines`, `jitter_rows_in_order`.
- `cobalt_files` (new): CobaltGcMedianFile.read (GcMedianReadDepth). Test `gc_median_lines`.
- `esvee_files` (new): EsveeDiscordantStats.read, DiscordantFragType. Test `disc_stats_lines`.
- `purple_files`: PurplePurity (every column, enum names checked as valueOf), `read_purple_purity`,
  `purple_purity_from_lines`, `purple_purity_filename`. Test `reads_every_column`.
- `lib.rs`: the four modules registered and listed.

## Open points

- The report's PNG text and clipping are hmf-plot-rs's (whole-pixel glyphs and clips; see "The report, drawn
  natively"); matching Cairo there needs sub-pixel glyph positioning and partial-coverage clipping in the
  rasteriser.

- SBX and ULTIMA: only the BQR quality bins depend on the sequencing type; checked with HG008's (Illumina)
  BQR file run as SBX and ULTIMA and with synthetic BQR rows (DUAL / SINGLE read types, MEDIUM qualities),
  not with real SBX or Ultima REDUX output (none at hand).
- Math.log10 on x86-64 (see above): the Java reference here is aarch64; a Linux x86-64 Java run could differ
  in the last bit of an MS indel score, which the 8-decimal output would show only rarely.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 11 on Google Batch ran qsee as the port, in Java mode, with every other step from the cache. The
image has no R, so its wrapper added `-skip_plot`. The three data files (gene limit of detection, status, vis
data) are identical to the first run's (qsee 1.1-beta.1); the R report was not drawn. ORANGE ran on them and
reported no qsee plot. qsee took 2.6 s and 21 MB, against Java's 38 s and 1.1 GB (with R). Since 2026-11-27
qsee-rs draws the report; for the pipeline to get it the image's `qsee` wrapper must stop adding `-skip_plot`
(it added it whenever Rscript was missing), after which ORANGE gets its qsee page as with Java.
