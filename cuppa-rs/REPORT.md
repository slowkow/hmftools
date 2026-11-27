# CUPPA in Rust: record of the port

This report covers how `cuppa-rs` was checked against Java CUPPA 2.5 and pycuppa, the faults found, and the
timings. `cuppa-rs` ports both CuppaDataPrep (Java) and pycuppa's `cuppa.predict` (Python), so it needs neither.
Finding ids: CU. Options, outputs and tests are in [README.md](README.md).

## Python reference environment

A private virtual environment made with uv 0.11.15 from CPython 3.11.10, with pycuppa installed from a copy of
its sources: numpy 1.26.4 (its wheel's OpenBLAS 0.3.23.dev, core armv8), pandas 2.0.3, scikit-learn 1.3.0, scipy
1.17.1, joblib 1.6.0. pycuppa's plot step needs R's ggh4x, installed in a private R library, so that the
reference runs to the end (exit 0) and draws its PNG. References run with OPENBLAS_NUM_THREADS=1 (2 threads give
the same bytes: OpenBLAS splits dgemv by output on macOS).

## Java reference

This repository's `cuppa/target/cuppa-2.5-jar-with-dependencies.jar` (sources following upstream ca04e6c39d of
2026-11-20). Since
upstream 3271b2ff48 (2026-11-17) it reads Viridian's report (`-viridian_dir`) and rejects the pipeline's
`-virus_dir` ("unregistered config item"). For the pipeline's command the reference is a
harness: the parent commit's DriverPrep and PrepConfig (Virus Interpreter's annotated file, the
2.5.1 image's behaviour; the only CUPPA sources that differ) compiled against the jar and put before it on the
class path. On HG008 the harness writes the pipeline's `cuppa_data.tsv.gz` byte for byte (gzip included), and on
the private pair that run's file. Java ran with `-Xmx8g -XX:ActiveProcessorCount=2`.

## How it was checked

A regression script (kept outside this repository with the outputs) runs, for each case, the Java prep (jar or
harness), then `python -m cuppa.predict`, then cuppa-rs with `-java_compat` (and again in the default mode), one
run at a time. It compares every output file (gzip files after decompression, and their bytes), the set of files
and the exit statuses, with `/usr/bin/time -l` per run; timings are medians of 3.

Inputs:
- GIAB HG008 (public): the pipeline's PURPLE, LINX and Virus Interpreter outputs, with the pipeline's command.
- The private pair: renamed copies (TUMOR / NORMAL) of the all-Java run's files
  (the ids found from the file names and VCF header, never printed).
- Synthetic samples built from HG008's public files:
  SYN1 (male; ALB and SFTPB indels, a filtered indel, an SNV with N in its context, an MNV, a `C>*` record, a
  FILTER "." record; no LINX catalog; an override fusion CBFB_MYH11, PROMISCUOUS_5 HIGH and PROMISCUOUS_3 LOW
  fusions; viruses in both formats), SYN2 (MSI; duplicate AMP rows, PARTIAL_AMP, GERMLINE_AMP, GERMLINE_MUTATION,
  NOT_REPORTED and a 0.12344 likelihood), SYN3 (generic variants CSV), SYN4 (empty inputs), a 12-sample cohort.
- Synthetic RNA: Isofox gene data for the model's 5,927 genes and 200 others
  (random TPMs), alternative splice junctions at 4,000 random reference sites, 50 off the sites, a duplicate.
- A cross-validation predictions file written by pycuppa.

### Results (Java mode)

| Case | What | Java + pycuppa | cuppa-rs | Result |
|---|---|---|---|---|
| hg008_pipeline | pipeline commands | prep 0.91 s 0.223 GB; predict 7.01 s 0.472 GB | prep 0.09 s 0.031 GB; predict 0.13 s 0.117 GB | identical (3 files) |
| t0_pipeline | private pair, pipeline commands | 0.65 s 0.202 GB; 7.26 s 0.471 GB | 0.06 s 0.016 GB; 0.13 s 0.118 GB | identical |
| cohort2 | HG008-T + private tumour, wildcards, 2 threads | 0.98 s 0.381 GB; 10.11 s 0.344 GB | 0.22 s 0.019 GB; 0.13 s 0.12 GB | identical |
| cohort12 | 12 samples, 2 threads, `--force_plot` | 3.12 s 2.094 GB; 45.34 s 0.408 GB | 1.14 s 0.070 GB; 0.16 s 0.124 GB | identical |
| cohort12_by_category | `-write_by_category -output_id`, directory input, `--sample_id` of 12 | 3.23 s 2.096 GB; 7.67 s 0.481 GB | 1.25 s; 0.14 s | identical |
| rna_all | ALL categories with synthetic RNA, `--clf_group all` | 1.08 s 0.254 GB; 8.32 s 0.503 GB | 0.15 s 0.048 GB; 0.13 s 0.122 GB | identical (combined, RNA rows) |
| rna_group_rna / rna_group_dna | `--clf_group rna` (no vis data) / dna | 1.07 s; 2.58 s / 7.06 s | 0.26 s; 0.37 s | identical |
| rna_cohort | 3 samples: RNA, NO_RNA, RNA file missing | 1.36 s 0.612 GB; 14.48 s | 0.34 s 0.056 GB; 0.38 s | identical |
| syn1_virus / syn1_viridian | harness `-virus_dir` / jar `-viridian_dir` | 0.88 s; 7.74 s | 0.20 s; 0.38 s | identical |
| syn2_msi, syn4_empty | MSI and drivers; empty inputs | 0.71 / 0.42 s; 7.36 / 6.75 s | 0.12 / 0.02 s; 0.37 s | identical |
| syn3_generic_snv / _dna | generic variants file, SNV only / DNA (Java NPE, CU3) | 0.41 s / exit 1 | 0.03 s / exit 1 | identical |
| hg008_categories | `-categories snv;sv;driver` (no sex feature) | | | identical |
| hg008_v37 | `-ref_genome_version 37` on GRCh38 data | | | identical |
| hg008_sample_data_dir | `-sample_data_dir` only, `V38`, `-log_debug` (jar) | 0.89 s | 0.20 s | identical |
| hg008_rename_compress | `--rename_classes`, `--compress_tsv_files`, `--log_to_file` | | | identical content (pandas' gzip has a time stamp) |
| hg008_no_sample_id, hg008_bad_group | single sample without `--sample_id`; `--clf_group combined` | exit 1 | exit 1 | identical |
| cv_predictions | `--cv_predictions_path` for one of 2 samples | 10.28 s 0.354 GB | 0.12 s 0.115 GB | identical |
| cohort2_sample_id | `--sample_id` of a cohort, `--log_path`, `--log_format` | 7.10 s | 0.36 s | identical |
| runner_syn1 | the jar's main class with `-python_path` (Java runs Python) | 8.74 s 0.493 GB | 0.25 s 0.111 GB | identical |
| err_* (6) | no sample, bad category, bad version, missing purity file, unknown option, missing output dir | exit 1 | exit 1 | identical (no files) |

Times are single runs except the first four rows and rna_all (medians of 3). pycuppa's times include its R plot
(about 1-2 s of the 7 s; 5.5 s are unpickling the classifier).

### Default (corrected) mode

Identical to Java and pycuppa in 20 of the 30 cases. The other 10 differ exactly where a correction applies:
the cohorts cohort2, cohort2_sample_id, cohort12, cohort12_by_category (CU7; each sample's rows are then byte-identical to its single-sample prediction, checked on cohort2 and
on the 12 copies of HG008-T in cohort12, except the one with an override fusion), rna_all and rna_group_rna (CU8),
rna_cohort (CU6 and CU7),
syn2_msi (CU6: the GERMLINE_AMP row), syn3_generic_dna (CU3: features written), hg008_categories (CU2: unknown
sex). cv_predictions is identical: its one computed sample is predicted alone in both modes.

## Speed and memory

- Prediction: 0.13 s and 0.12 GB against pycuppa's 7 s and 0.47 GB. Most of it is unpacking the embedded model
  (10 MB of gzip, 21 MB, 2.4 million floats) and building its tables; the features file, the arithmetic and the
  outputs take a few ms. With the pickle as `--classifier_path` the port reads only its size and last 8 bytes.
- Feature extraction: 0.06-0.09 s and 16-31 MB against Java's 0.65-0.91 s and 0.2 GB. Of HG008's 0.09 s, 64 ms
  read the somatic VCF (3.6 MB of BGZF, 19,254 records; inflated with hmf-common's zlib-rs reader on a second
  thread), 11 ms the SV VCF and LINX files, the rest the drivers (the VCF is read once for SNV and DRIVER; Java
  reads it twice). Cohorts: 1.1 s for 12 samples at 2 threads against Java's 3.1 s and 2.1 GB.

## Model resource

`python -I cuppa-rs/scripts/export_model.py cuppa_classifier.38.pickle.gz cuppa-rs/resources/cuppa_model.38.bin.gz`
in the environment above (9.85 MB; deterministic: sha256 deea3b4d2625ccd75109ac0aeb9662870ec7662d1c162a1199bfb529281760f2
from the pickle with sha256 3e0f9490a9f37ac252b23b38b6b991245b72bcdf8f963b2a91ade82547c9fe70). It checks the
versions (scikit-learn 1.3.0, numpy 1.x) and the classifier's structure (step names, column selectors, options)
and writes every parameter the prediction uses: the required features, noise profiles, class profiles and their
L2 norms as numpy computes them at each prediction (sample-independent, so stored rather than reproduced: their
summation order depends on the array layout, and the float32 alt-SJ profiles are summed in float32), the logistic
regressions (coef, intercept, the event features' training means), max-scaler maxima, chi2-selected features,
isotonic calibrators (thresholds), fusion overrides, sex keywords, the signature quantile transformers, the
cross-validation performance, and the source pickle's SHA-256, size and last 8 bytes. Floats keep their IEEE
bytes (float32 stays float32), byte-shuffled for compression, gzip -9 with mtime 0. `--classifier_path` takes the
pickle it was exported from or an exported file; another pickle must be exported first.

## Numerics (how pycuppa's numbers are reproduced)

- pycuppa's last 2-3 of 17 digits depend on the machine: the pipeline's HG008 run (Linux x86-64) and the local
  reference (macOS arm64) differ there in 117 probabilities and 151 odds of vis_data (the private pair: the same).
  The port reproduces the macOS arm64 reference (numpy 1.26.4 wheels with OpenBLAS 0.3.23, armv8 kernels).
- Products: `x @ coef.T` and `np.dot(x, profiles)` go to OpenBLAS dgemv_t (kernel/arm64/gemv_t.S: 8 fused partial
  sums over blocks of 32, folded ((v1+v2)+v3)+v4, pairs of 4, a fused tail); with float32 coefficients numpy
  widens a row-major copy and takes dgemv_n (one fused chain); for several samples dgemm (K blocks of 128, the last
  two balanced and rounded up to 8, each a fused chain, added to C). Checked on 1,306 random products, the fixture,
  and the cohort cases.
- Sums: numpy's pairwise summation for one row; pandas sums the rows of a cohort's frames column after column
  (determined on the cohort runs: cosine norms, calibration, sex filter). np.interp without fused multiply-add
  (checked on 20,000 values); exp, log1p, pow, log10 from the C library, as numpy's arm64 loops.
- Ties: np.argsort's introsort for doubles (the top classes and features); pandas' rank(method="first").
- Text: Python's float repr (numpy's str of a float64), pandas' CSV writing, Java's String.format for the features.

## Findings (CU)

Corrected by default, Java's or pycuppa's behaviour with `-java_compat`:
- **CU2 (pycuppa SexProbFilter): a missing sex feature is taken as male.** fill_missing sets a missing
  `event.trait.is_male` to -1e-8 and `astype(bool)` makes that True, so the female classes (ovary, endometrium,
  ...) are zeroed for a sample without SAMPLE_TRAIT features (`-categories` without it, or a cohort row missing
  it); the code's own rule is that an unknown sex leaves the probabilities as they are. Corrected: a sex feature
  missing (or NA) in the features file is unknown. Case hg008_categories.
- **CU3 (Java SomaticVariantsLoader.loadFromGenericFile): the type filter is inverted.**
  `if(variantTypes == null && !variantTypes.contains(...))` never filters by type and throws a
  NullPointerException when no types are given (DriverPrep's indel drivers of a microsatellite-stable sample read
  with `-somatic_variants_dir`): exit 1, no features. Corrected: the types filter as for the VCF, none given means
  all. Case syn3_generic_dna.
- **CU4 (Java CuppaDataPrep, single sample): missing RNA data crashes.** GeneExpressionPrep and
  AltSpliceJunctionPrep return null for a missing Isofox file (meant for cohorts, where the sample's column is
  left null), but the single-sample path calls `dataItems.addAll(null)`: NullPointerException, no features file.
  Corrected: that category is left out with a warning. Unit test cu4.
- **CU5 (Java PythonInterpreter.checkValid): only an interpreter named exactly "python" is accepted.**
  `path.getFileName().startsWith("python")` is Path.startsWith, which compares whole names, so `/usr/bin/python3`
  is an "Invalid python interpreter". The port needs no interpreter: `-python_path` is optional and unused by
  default; with `-java_compat` it is required and checked as Java does.
- **CU6 (hmf-common DriverType.isGermline): GERMLINE_AMP is not germline.** A GERMLINE_AMP row in the driver
  catalog becomes a somatic `<gene>.mut` feature (".mut" as it is not AMP). Not met in the pipeline's catalogs.
  Corrected: skipped like the other germline types. Cases syn2_msi, rna_cohort.
- **CU7 (pycuppa with numpy, OpenBLAS and pandas): a sample's numbers depend on the cohort it is predicted in.**
  For one sample numpy multiplies with dgemv and pandas sums rows pairwise; for several, dgemm and sequential row
  sums, and the fusion overrider renormalises every row once any sample has an override fusion: HG008-T's top
  dna_combined probability is 0.39593293399613794 alone and 0.3959329339961264 with one other sample. Corrected:
  every sample computed as if alone (its rows identical to its single-sample prediction).
- **CU8 (pycuppa / CuppaDataPrep): the alt-SJ classifier gets raw counts but was trained on ln(count + 1).**
  CuppaDataPrep writes the fragment counts; pycuppa's AltSjClassifier comments out its Log1pTransformer ("log
  transform in DataLoader beforehand"), and no step of the prediction takes the log. The model's maxima prove the
  training scale: the smallest is 0.693 = ln 2 (a count of 1), the largest 9.2 (about 10,000). Raw counts of 2 and
  more therefore saturate the max scaling at 1, and the alt-SJ classifier sees nearly binary features. RNA runs
  only (the pipeline runs DNA). Corrected: the counts are log-transformed when the features are read (the features
  file stays as Java writes it). Cases rna_all, rna_group_rna, unit test cu8.
- **CU9 (plot_predictions.R): a signature quantile of -Inf is drawn as "0", with an "NA" legend key.** pycuppa's
  SigCohortQuantileTransformer gives -(first quantile / value) below the cohort's range, which is -Inf when the
  sample has none of a signature that every sample of the class has (a sample with few SNVs). The script bins the
  quantiles with `cut(breaks = c(-Inf, -0.004, ...))`, whose first interval is open at -Inf, so -Inf becomes NA:
  the cell takes `na.value` (grey95, the colour of "0", in the expected range) and the legend gains a key "NA",
  while its label says "-Inf". Corrected: -Inf is in the first class, "<0 (below expected range)" (light cyan),
  and the legend has no NA key; the cell's label stays "-Inf". Case syn4_empty (empty inputs).

Not reproduced:
- **CU1 (Java PredictionRunner): the Python command is one string for `bash -c`.** The arguments are joined with
  spaces unquoted (ShellCommand), so a path or sample id with a space breaks the command and one with shell
  syntax is executed by bash. The port runs the prediction in process with the arguments as given, in both modes.

Kept as Java and pycuppa do, recorded:
- The jar rejects the pipeline's `-virus_dir`; the port takes it (Virus Interpreter's file, likelihood 0.5 for a
  LOW driver likelihood, as 2.5.1) as well as `-viridian_dir`. With neither, the jar's choice (Viridian's report
  from the sample data directory).
- The model's 115 chrY gen_pos features are never written by the 2.5 prep (pycuppa's loader rejects chrY) and are
  filled with -1e-8: harmless, as neither the class profiles nor the noise profile have mass on chrY.
- htsjdk types a `C>*` record as an SNP (its alt becomes the ref, it counts in gen_pos, and logs an invalid
  bucket for the 96 contexts); kept.
- DEL and other non-amplification drivers become `<gene>.mut`, as the model was trained.
- pycuppa logs some errors with arguments its format does not use (duplicated features, missing snv_count), and
  raises a TypeError instead of its CalledProcessError when R is missing; cosmetic, the port logs plain messages.

## The summary plot

`<sample>.cuppa.vis.png` (one sample) and `cuppa.vis.pdf` (several; one page per sample) are drawn natively
(`src/plot.rs`, with hmf-plot-rs's `ggp`, `gtable`, `patchwork`, `rdev` and `pdf_figure`), as pycuppa's
`plot_predictions.R` draws them in the pipeline's CUPPA 2.5.1 image (R 4.2.3, ggplot2 3.5.1, ggh4x 0.2.8, patchwork
1.2.0; Cairo png() with DejaVu Sans; that image's render of HG008's vis data is byte-identical to the pipeline's
PNG). `cuppa-rs predict` draws it where pycuppa does: not in RNA-only mode, a PDF for several samples, nothing
above 10 samples without `--force_plot`; `cuppa-rs plot <vis_data> <out>` draws it from a vis data file.

How it was checked (2026-11-27, the ports in Java mode): R drew each case in the image (docker, amd64 emulated),
and again with every grid viewport's device box dumped; the port wrote its grobs' boxes. Layout: the boxes paired
by name (gtable and grid's own levels dropped), then every box without a name partner looked for by position
among all of the other side's; the PNG compared pixel by pixel; the PDFs by pages and page size, `pdftotext
-layout`, every word's box (`pdftotext -bbox`) against R's, and the pages at 100 dpi.

| case | input | R in the image: s, MB | port: s, MB | PNG pixels differing by >10% / >25% / >50% | PDF |
|---|---|---|---|---|---|
| hg008 | HG008's vis data (DNA) | 8.20, 284 | 0.26, 221 | 108,755 / 50,249 / 7,030 of 19.8 M | |
| t0 | the private pair | 7.57, 278 | 0.24, 207 | 98,251 / 39,371 / 7,510 | |
| rna_all | DNA + synthetic RNA | 8.62, 292 | 0.30, 249 | 131,308 / 63,892 / 7,576 | |
| syn4_empty | empty inputs (CU9 in Java mode) | 7.29, 277 | 0.24, 206 | 105,403 / 46,531 / 11,335 | |
| syn2_msi | MSI and drivers | 7.86, 292 | 0.29, 240 | 109,543 / 47,871 / 5,953 | |
| hg008_categories | snv, sv, driver only | 7.32, 283 | 0.26, 220 | 109,058 / 50,313 / 7,030 | |
| syn1_virus | a virus | 8.38, 297 | 0.29, 250 | 118,727 / 50,754 / 7,069 | |
| rna_only | RNA only | stops (no `snv_count`) | stops, same message | | |
| cohort2 | HG008-T and a private tumour | 11.34, 195 | 0.03, 9 | | 2 pages, text identical, 2,267 words at 0.000 pt, 0 / 0 pixels |
| rna_cohort | DNA + RNA cohort | 17.27, 199 | 0.04, 11 | | 3 pages, text identical, 3,668 words, 0 / 0 / 0 pixels |
| cohort12 | 12 samples | 58.97, 207 | 0.13, 30 | | 12 pages, text identical, 14,053 words, 0 pixels on every page |

- Size: R's in every case (6000 x 3075 to 6000 x 3750 at 300 dpi; PDF pages 1440 x 900 pt).
- Layout (hg008, t0, rna_all, syn4_empty): the 185 (rna_all 192) elements paired by name are within 0.0000 px of
  R's. Of the rest, the port's 76 (133) panel boxes of the probability heatmap are cells of R's grid (R's dump has
  no viewports below that panel group, but its rows and columns are R's axes and strips), and R's unpaired boxes
  are 248 (214) of zero size, 251 (259) axis viewports sharing three edges with the port's cell (an axis's own
  viewport, as wide as its labels), and 26 (28) label columns inside axes, which the port does not record.
- Visible differences (PNG; none at normal size): text is rasterised by hmf-plot-rs with whole-pixel glyph
  positions (Cairo positions glyphs to a fraction of a pixel; the bold titles are BUI's synthetic bold); panels
  clip to whole pixels, so a border drawn inside a clipped panel differs in its outermost pixel line (9% of the
  pixels differing by >25% on hg008, none by >50%); anti-aliasing of lines and tiles is BUI's. Nothing is drawn
  in another place, colour or size.
- PDF: identical to R's at 100 dpi on every page (after three fixes on 2026-11-27: legend key backgrounds,
  rectangle outlines written as R writes them, and PDF dash lengths; hmf-plot-rs README).
- Determinism: a second run gives the same bytes; built on Linux arm64 (`rust:1.97-bookworm`) the HG008 PNG has
  the same bytes as on macOS arm64.
- Speed and memory against native R on the same Mac (one thread, three runs each): the HG008 PNG 0.26 s and
  220 MB against 4.3-4.8 s and 432 MB; the 12-page cohort PDF 0.12 s and 29 MB against 41 s and 400 MB. The
  port's memory is mostly the 6000 x 3300 canvas.
- ORANGE: with the port's plot in `cuppa/`, ORANGE's JSON has `cuppaSummaryPlot` as Java's (the whole JSON
  byte-identical to Java's on HG008) and its CUPPA chapter shows the plot; without it (before this) the JSON had
  null there and the chapter had no plot.

## Open points

- **The platform's exp:** the softmax takes `f64::exp` from the platform's math library, as numpy does, so the
  port matches pycuppa on macOS byte for byte; built on Linux (glibc) its last digits can differ, as pycuppa's own
  do between machines (the HG008 test compares within a relative 1e-12 there). Inside the Linux pipeline, the
  pipeline's pycuppa also used x86-64 BLAS kernels, whose summation order differs again.
- **rustc 1.100.0:** a release build with rustc 1.100.0 (Linux, arm64 and x86-64) typed every SNV of the Java
  tests' VCF as UNDEFINED in `hts_type` (plain safe code; the same build with a debug print in that function
  typed them correctly), so four prep tests failed. With rustc 1.97.1, the version all checks used, they pass on
  Linux and macOS. The image is pinned to rustc 1.97 until the cause (a compiler fault, or undefined behaviour
  elsewhere in the process) is found.

- **Other machines:** the port's numbers are those of pycuppa on macOS arm64. On Linux x86-64 (the pipeline)
  pycuppa's last digits differ (OpenBLAS's Haswell/SkylakeX kernels, numpy's AVX512 exp/log/pow); the port gives
  its macOS numbers there too, apart from the C library's exp, log1p, pow and log10 (glibc and Apple's libm agree
  in all but rare last-bit cases; not checked on Linux).
- **Plots:** drawn natively since 2026-11-27 (see "The summary plot" above); text in the PNG is rasterised by
  hmf-plot-rs (whole-pixel glyphs), and panels clip to whole pixels, so pixels at glyph and panel edges differ
  from R's.
- **Other classifiers:** need the exporter (Python with scikit-learn 1.3.0) once; a pickle reader in Rust (as
  chord-rs reads .rds) would remove that.
- **hmf-common-rs:** `variant_context::is_symbolic_allele` (and so `VcfRecord::variant_type`) calls `*`
  symbolic, where htsjdk 2.23 types it as a base (SNP or INDEL); CUPPA types records itself. Other ports reading
  `*` alleles may differ from htsjdk.
- **oncoanalyser:** the CUPPA module's versions.yml runs `java`, `R`, `Rscript` and `python` to report versions;
  running the port needs that part changed, as for the other ports.

## Shared-code additions (hmf-common-rs)

- `sigs` (new module): LeastSquaresFit (`sigs/lsq.rs`) and Matrix / loadMatrixDataFile (`sigs/matrix.rs`) moved
  from sigs-rs unchanged (crate paths only); SnvSigUtils (bucket names, variantContext, swapDnaBase, the COSMIC
  signatures resource `sigs/snv_cosmic_signatures.csv`, a copy of hmf-common's), PositionFrequencies' static
  functions (standard lengths, initialisePositionCache with Java's start rule as a parameter, getBucketIndex,
  getChromosomeFromIndex, getPositionFromIndex) and `PositionCounts` (buildMap false). Tests in the module.
  sigs-rs could use it in place of its copies.
- `config_builder::convert_wildcard_sample_path` and `java_replacement` (ConfigUtils.convertWildcardSamplePath with
  Matcher's replacement rules), from chord-rs's copy.

Used from hmf-common-rs: `config_builder`, `compat`, `log`, `fmt` (Java's %.nf / %e), `java_text`, `java_hash`
(the reportable signatures' HashMap order), `variant_context` (record parsing), `sv` (StructuralVariantFactory,
EnrichedStructuralVariantFactory, StructuralVariantData), `purple_files` (driver catalog), `gz` (zlib-rs
inflating), `regions`, `ref_coords`, `sigs`.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 11 on Google Batch ran CUPPA as the port (the image's `cuppa`, and its `python` stand-in for
`python -m cuppa.predict`), in Java mode, with every other step from the cache. The features file is identical
to the first run's. The prediction summary and vis data differ from the pipeline's pycuppa (Linux x86-64) only in
the last digits: 10 and 247 values, at most 2.4e-13 relative. ORANGE, which ran again on them, carried those
digits into its JSON and reported no CUPPA plot (the port drew none then; it does since 2026-11-27, see "The
summary plot"); everything else in the JSON was unchanged.
CUPPA took 1.7 s and 21 MB, against Java and pycuppa's 16 s and 494 MB. The rerun's tasks cost $0.02.
