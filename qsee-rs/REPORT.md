# qsee in Rust: record of the port

This report covers how `qsee-rs` was checked against Java qsee 1.0, the faults found in the Java code, and the
timings. `qsee-rs` is the QSEE_METRICS step of oncoanalyser: QC features from the PURPLE, BamTools, COBALT, ESVEE
and REDUX outputs, QC statuses against thresholds, the gene limit-of-detection file, and the R report. The Java
reference is this repository's `qsee/target/qsee-1.0-jar-with-dependencies.jar` (sources following upstream
ca04e6c39d of 2026-11-20), on Java 23. `resources/plot_qc.R` is the jar's R script, unchanged. Options and
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
| lod: LimitOfDetectionCalcs, GeneLodCalculator, GeneLodFile, VariantDetectabilityStatus | `lod.rs` |
| cohort: FeatureMatrix, PercentileTransformer, FeaturePercentiles, CohortPercentiles(File), CohortPercentilesTrainer, CohortFeaturesWriter | `cohort.rs` |
| plot: VisDataFile, QseePlot, QseePlotTask (RExecutor) | `vis.rs` |
| hmf-common: DelimFileReader/Writer, createBufferedWriter (gzip), BqrFile, DuplicateFrequency, JitterCountsTableFile, CobaltGcMedianFile, EsveeDiscordantStats, PurplePurity (full), PurityContextFile.readWithQC | hmf-common-rs additions (below) and `prep.rs` |

Not ported: PDFBox's merging of the per-sample PDFs for `-merge_plots` (multi-sample mode only; the port
keeps the per-sample PDFs and warns), and `compare_cohorts.R`, which no Java class runs.

### The plot step

The port runs the jar's `plot_qc.R` with the 8 arguments QseePlotTask passes (`Rscript`, the script from a
temporary file, as RExecutor), unless `-skip_plot`. The R step reads only the vis data file and the cohort
percentiles file. On the two pipeline cases the reports drawn from the port's files equal Java's: the PNG
byte for byte, the PDF apart from its CreationDate and ModDate. Panels and their data: `qsee-rs/README.md`.
R here: `/usr/local/bin/Rscript` with dplyr, ggplot2, patchwork and scales already installed (nothing was
installed; the scripts set R_PROFILE_USER and R_ENVIRON_USER to /dev/null).

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

Remarks, not findings:
- The discordant fragment features divide fragment counts by total reads (ESVEE notes its counts are in
  fragment terms, "hence see doubling in rate calc"), so they are half the share of reads; the cohort
  percentiles use the same definition, so the comparison is consistent and nothing is changed.
- `-output_id` is registered with the description "Output directory" (QseeConfig), as `-help` shows.
- `-num_percentiles 1` gives one percentile NaN (100 / 0 steps); the file has a column `Pct_NaN`.
- `compare_cohorts.R` is in the jar but no Java class runs it.

## Speed and memory

The port is 14 (HG008) and 16 (private pair) times faster than Java on the pipeline's data files, with a
fortieth of the memory: 0.21 s and 34 MB against 2.92 s and 1,359 MB (HG008), 0.27 s and 38 MB against 4.22 s
and 1,433 MB (private pair), both with `-skip_plot`. With the R report (the pipeline's command) both take
about 9 s more, R's time: 9.0 s against 12.3 s (HG008), 9.6 s against 14.2 s; R's peak (about 330 MB) is then
the port's. In multi-sample mode (12 synthetic samples) 1.0 s against 13.8 s.

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

- `-merge_plots` (multi-sample PDF merging with PDFBox) is not ported; the per-sample PDFs are kept.
- SBX and ULTIMA: only the BQR quality bins depend on the sequencing type; checked with HG008's (Illumina)
  BQR file run as SBX and ULTIMA and with synthetic BQR rows (DUAL / SINGLE read types, MEDIUM qualities),
  not with real SBX or Ultima REDUX output (none at hand).
- Math.log10 on x86-64 (see above): the Java reference here is aarch64; a Linux x86-64 Java run could differ
  in the last bit of an MS indel score, which the 8-decimal output would show only rarely.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 11 on Google Batch ran qsee as the port, in Java mode, with every other step from the cache. The
image has no R, so its wrapper added `-skip_plot`. The three data files (gene limit of detection, status, vis
data) are identical to the first run's (qsee 1.1-beta.1); the R report was not drawn. ORANGE ran on them and
reported no qsee plot. qsee took 2.6 s and 21 MB, against Java's 38 s and 1.1 GB (with R).
