# cuppa-rs

A Rust port of hmftools CUPPA 2.5 (`../cuppa`), cancer of unknown primary prediction, in both of its parts:
the Java feature extraction (`com.hartwig.hmftools.cup.prep.CuppaDataPrep`) and pycuppa's prediction
(`python -m cuppa.predict`, scikit-learn 1.3.0), so that neither Java nor Python is needed. The classifier is
exported once from its pickle (`scripts/export_model.py`) into `resources/cuppa_model.38.bin.gz`, which the
binary embeds. With `-java_compat` the output files are byte-identical to Java's and pycuppa's; faults found in
both are corrected by default (findings CU1-CU9 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/cuppa-rs
```

The two commands of oncoanalyser's CUPPA task (commands/CUPPA.command.sh), with the binary in place of the
`cuppa` wrapper and of Python:

```
cuppa-rs -Xmx12g com.hartwig.hmftools.cup.prep.CuppaDataPrep -sample TUMOR -categories DNA \
    -purple_dir purple -linx_dir linx_somatic -virus_dir virusinterpreter -ref_genome_version 38 -output_dir cuppa/
cuppa-rs predict --sample_id TUMOR --features_path cuppa/TUMOR.cuppa_data.tsv.gz --clf_group DNA \
    --classifier_path cuppa_classifier.38.pickle.gz --output_dir cuppa/
```

- `cuppa-rs prep` is the same as the class name; `cuppa-rs -m cuppa.predict` the same as `predict`.
- Without a command (or with `com.hartwig.hmftools.cup.cli.PredictionRunner`) it is the jar's main class: the
  features (unless `-features_path` is given) and then the prediction in the same process, with the arguments
  Java passes to Python.
- JVM options (`-Xmx...`, `-XX:...`, `-D...`) are accepted and ignored; `cuppa-rs -version` prints
  `Cuppa version 2.5`.
- `cuppa-rs plot <sample>.cuppa.vis_data.tsv <out.png|out.pdf>` draws the summary plot from a vis data file
  (as `Rscript plot_predictions.R`).

## Options

**Feature extraction** (as the jar's PrepConfig): `-sample` (and `-rna_sample`) or `-sample_id_file` (TSV with
`SampleId` and optionally `RnaSampleId`; blank means the SampleId, `NO_RNA` skips the RNA categories);
`-categories` (DNA, the default; RNA; ALL; or a `;` list of SNV, SV, SAMPLE_TRAIT, DRIVER, GENE_EXP, ALT_SJ);
directories `-sample_data_dir`, `-purple_dir`, `-linx_dir`, `-viridian_dir`, `-isofox_dir`,
`-somatic_variants_dir` (generic `<sample>.somatic_variants.csv` in place of PURPLE's VCF), with `*` replaced by
the sample id; `-virus_dir` (Virus Interpreter's `<sample>.virus.annotated.tsv`, as CUPPA 2.5.1 and the pipeline;
the jar itself reads Viridian's `<sample>.viridian.report.tsv` from `-viridian_dir` or the sample data directory);
`-ref_genome_version` (37 or 38, V37 by default), `-ref_alt_sj_sites` (needed for ALT_SJ), `-output_dir`,
`-output_id`, `-write_by_category`, `-threads` (samples at a time in a cohort), `-log_level`, `-log_debug`.
The PredictionRunner adds `-features_path`, `-classifier_path`, `-clf_group`, `-cv_predictions_path` and
`-python_path` (not needed; checked as Java checks it only with `-java_compat`, CU5).

**Prediction** (as pycuppa's argparse options): `--classifier_path` (the pickle the embedded model was exported
from, recognised by its size and gzip trailer, or a file written by `scripts/export_model.py`),
`--features_path` (a single- or multi-sample features file, or a directory of `cuppa_data.cohort.<category>.*`
files), `--output_dir`, `--sample_id`, `--clf_group` (all, dna or rna), `--cv_predictions_path` (samples found in
a CuppaPrediction file take its predictions), `--compress_tsv_files`, `--force_plot`, `--rename_classes
A=B,C=D`, `--log_to_file`, `--log_path`, `--log_format`.

**Java and pycuppa behaviour:** `-java_compat` (anywhere on the command line, or `CUPPA_JAVA_COMPAT=1`) for every
finding; `-java_compat_except CU2,...` (or `CUPPA_JAVA_COMPAT_EXCEPT`) for all but those listed.

## Outputs

- `<sample>.cuppa_data.tsv.gz` (or `cuppa_data.cohort[.<category>][.<output_id>].tsv.gz` for several samples):
  Source, Category, Key, Value per feature, gzip as java.util.zip writes it (byte-identical);
- `<sample>.cuppa.pred_summ.tsv`: per classifier the top 3 classes and probabilities, with the top signatures
  (snv96) and the top event features per class (event), as pandas writes them;
- `<sample>.cuppa.vis_data.tsv` (not in RNA mode): probabilities, signature quantiles, event feature odds and the
  cross-validation performance, with ranks;
- `<sample>.cuppa.vis.png` (or `cuppa.vis.pdf` for a cohort, one page per sample): the summary plot, drawn
  natively (`src/plot.rs`, with hmf-plot-rs) as pycuppa's R script `cuppa/visualization/plot_predictions.R`
  draws it in the CUPPA 2.5.1 image (R 4.2.3, ggplot2 3.5.1, ggh4x, patchwork), at R's pixel size. Every element
  is where R puts it (within 0.01 px of R's viewports); text is DejaVu Sans as there, but rasterised by
  hmf-plot-rs (whole-pixel glyph positions), so pixels at glyph edges, panel edges (whole-pixel clipping) and
  anti-aliased lines differ from Cairo's. The PDF has R's pages, page size, text and word positions (Helvetica,
  R's kerning) and draws as R's at 100 dpi to the pixel. `cuppa-rs predict` draws it where pycuppa would (not in
  RNA-only mode; above 10 samples only with `--force_plot`). R is not needed. Checks: [REPORT.md](REPORT.md),
  "The summary plot".

## How the prediction is reproduced

The prediction follows each scikit-learn and pycuppa step (noise profiles, cosine similarities to the class
profiles, non-best scaling, max scaling, logistic regressions, isotonic calibration, fusion override, sex filter,
the DNA x RNA combination, signature quantiles, feature contributions) with numpy's arithmetic: OpenBLAS 0.3.23's
armv8 kernels for the products (dgemv_t, dgemv_n for float32 coefficients, dgemm for cohorts), numpy's pairwise
summation, np.interp, np.argsort's introsort for ties, pandas' row sums, Python's float repr. Its numbers equal
pycuppa's on macOS arm64 (numpy 1.26.4 wheels); pycuppa itself gives other last digits on Linux x86-64 (other
kernels and SIMD math). By default each sample of a cohort is computed as if alone (CU7).

## Tests

`cargo test --release -- --test-threads=2`: the Java unit tests of the prep (`src/tests/java_*`, with the Java
test resources in `testdata/java`), pycuppa's tests of the prediction components (`src/tests/pycuppa.rs`), and
the port's tests (`src/tests/port.rs`: numpy's arithmetic against `testdata/numpy_fixture.tsv`, written by
numpy through `scripts/numpy_fixture.py`; the embedded model; HG008's prediction against pycuppa's in both modes;
each finding in both modes; `src/plot.rs`: stringi's word wrapping and R's number formats as the plot uses them).
The plot's layout model is tested in hmf-plot-rs (`ggp`, `patchwork`, `pdf_figure`). The comparisons with Java,
pycuppa and R are in [REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
