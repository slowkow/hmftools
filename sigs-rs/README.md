# sigs-rs

A Rust port of hmftools SIGS 1.3 (`../sigs`): a sample's passing SNVs counted in the 96 trinucleotide
buckets and fitted to SNV signature definitions (by default the COSMIC 30) by non-negative least squares,
with the counts the fit does not explain reported as `MISALLOC`. It takes the same options and writes the
same files as the Java `SignatureFitter`; with `-java_compat` they are byte-identical to Java's. Faults found
in the Java code are corrected by default (findings SI1-SI7 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/sigs-rs
sigs-rs -sample TUMOR -somatic_vcf_file TUMOR.purple.somatic.vcf.gz -signatures_file snv_cosmic_signatures.csv \
    -log_level INFO -output_dir sigs/
```

This is the command line oncoanalyser runs (SIGS). JVM options (`-Xmx...`, `-XX:...`) are accepted and
ignored, so the bioconda wrapper's command line works, and `sigs-rs -version` prints `Sigs version 1.3`.

## Options

As Java's (`sigs-rs -help` lists them):
- **Samples:** `-sample`, or `-sample_id_file` (a CSV with a `SampleId` column, or one column; lines starting
  with `#` are skipped); more than one sample writes one fit file for all.
- **Inputs:** `-somatic_vcf_file` (bgzipped or plain, indexed or not), or `-sample_counts_file` (a CSV of
  bucket counts: `BucketName` and a column per sample); `-signatures_file` (CSV, a column per signature, 96
  rows; without it, the COSMIC signatures hmf-common carries, in `resources/`).
- **Filters** (each applied only when given): `-subclonal_min`, `-subclonal_max` (SUBCL), `-ploidy_min`,
  `-ploidy_max` (PURPLE_VCN), `-tumor_only` (drops LOW somatic likelihood, hotspots, a non-zero GND_FREQ
  and allele frequencies below 5%).
- **Fit:** `-min_alloc` (0.01) and `-min_alloc_perc` (0.0005): smaller allocations are not written.
- **Position frequencies:** `-position_bucket_size` (0: none), `-write_pos_freq_coords`,
  `-ref_genome_version` (37 by default, as in Java, or 38).
- **Output:** `-output_dir`, `-output_id`; `-log_level`, `-log_debug`.
- **Java behaviour:** `-java_compat` gives Java's behaviour for every finding; `-java_compat_except SI1,...`
  for all but those listed. The environment variables `SIGS_JAVA_COMPAT=1` and `SIGS_JAVA_COMPAT_EXCEPT` do
  the same.

## Outputs

- `<sample>.sig.snv_counts.csv`: the 96 bucket counts (one sample);
- `<sample>.sig.allocation.tsv`: signature, allocation (`%.3f`) and percent (`%.5f`) of each signature
  above both minimums, then `MISALLOC` (one sample);
- `sig_snv_fit.csv`: `SampleId,Signature,Allocation,Percent` (Java's `%g`) for more than one sample;
- `<sample>.sig.pos_freq_counts.csv` with `-position_bucket_size`.

With `-output_id X` the files are `<sample>.X.sig.allocation.tsv`, `<sample>.sig.snv_counts.X.csv` and
`sig_snv_fit.X.csv`, as in Java.

## How it works

The VCF is read once; filtered records and contigs other than 1-22, X and Y are skipped before a record is
parsed. For the others, what SmallVariantFactory and SIGS read is taken from the record (allelic depths,
TNC, SUBCL, PURPLE_VCN, SOM_LH, HOTSPOT, GND_FREQ). The fit is LeastSquaresFit (Lawson and Hanson's NNLS)
step for step, so the contributions are Java's to the last bit; rounding uses StrictMath's `pow(10, -n)`,
which for 10^-5 is a unit in the last place below the literal `1e-5`, as Java's `Math.pow` is.

## Tests

`cargo test --release` runs the Java unit test of SIGS (SignaturesTest, `src/tests/signatures_test.rs`) and
tests of the port (`src/tests/port.rs`: buckets and contexts, the fit, rounding, the file formats, the matrix
reader, and the corrected findings in both modes). The comparisons with Java are in [REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
