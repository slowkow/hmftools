# chord-rs

A Rust port of hmftools CHORD 2.1.2 (`../chord`): homologous recombination deficiency (HRD) prediction
from a tumor's mutation contexts. Java CHORD counts the contexts and then runs an embedded R script
(`chord_predict.R`) that predicts with a randomForest model (`CHORD.rds`); the port does both in Rust and
needs no R. It takes the same options and writes the same files as the Java `ChordApplication`; with
`-java_compat` the files are byte-identical to Java's. Faults found in the Java code and the R script are
corrected by default (findings CH1-CH10 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/chord-rs
chord-rs -sample TUMOR -snv_indel_vcf_file TUMOR.purple.somatic.vcf.gz -sv_vcf_file TUMOR.purple.sv.vcf.gz \
    -ref_genome ref.fasta -output_dir chord/
```

These are the options oncoanalyser passes (commands/CHORD.command.sh). JVM options (`-Xmx...`, `-XX:...`)
are accepted and ignored, so the bioconda wrapper's command line works, and `chord-rs -version` prints
`CHORD version 2.1.2` as oncoanalyser expects.

## Options

As Java's (`chord-rs -help` lists them):
- **Samples:** `-sample` (one id, or several separated by ','), or `-sample_id_file` (a CSV with a
  `SampleId` column, or one column; lines starting with '#' are skipped). Several samples make a cohort
  run, written to `cohort.chord.*.tsv`, with `-threads` samples at a time.
- **Inputs:** `-snv_indel_vcf_file` and `-sv_vcf_file` (PURPLE's somatic and SV VCFs; `*` in a path is
  replaced by the sample id), or `-purple_dir` (`<dir>/<sample>.purple.somatic.vcf.gz` and
  `.purple.sv.vcf.gz`); `-ref_genome` (indexed FASTA).
- **Counting:** `-include_non_pass` (count filtered variants too), `-write_detailed_files` (a line per SNV,
  INDEL and SV with its context: `<sample>.chord.{snv,indel,sv}.details.tsv`).
- **Model:** `-chord_model_file`: a randomForest model saved by R (`.rds`, read without R) or a forest
  exported by `scripts/export_forest.R`; without it, the embedded CHORD model.
- **Output:** `-output_dir`, `-output_id` (a suffix before `.tsv`), `-log_level`, `-log_debug`.
- **Java behaviour:** `-java_compat` gives Java's behaviour for every finding; `-java_compat_except CH1,...`
  for all but those listed. The environment variables `CHORD_JAVA_COMPAT=1` and `CHORD_JAVA_COMPAT_EXCEPT`
  do the same.

The other two Java entry points are commands: `chord-rs prep <options>` (ChordDataPrep: the counts
only) and `chord-rs predict -mut_contexts_file <tsv> -output_file <tsv> [-chord_model_file <model>]`
(ChordModel: the prediction only). The Java class names work in place of `prep` and `predict`.

## Outputs

- `<sample>.chord.mutation_contexts[.<id>].tsv`: the sample's counts of 96 SNV contexts (A[C>A]A, ...),
  30 INDEL contexts (del.rep.len.1, ..., ins.none.len.5) and 19 SV contexts (DEL_0e00_1e03_bp, ..., TRA);
- `<sample>.chord.prediction[.<id>].tsv`: p_BRCA1, p_BRCA2, p_hrd, hr_status (HR_deficient,
  HR_proficient or cannot_be_determined), hrd_type (BRCA1_type, BRCA2_type or none) and the QC remarks
  (MSI, fewer than 100 INDELs, fewer than 30 SVs), as R's write.table writes them;
- with `-write_detailed_files`, the details files above.

## How it works

The somatic VCF is read once (Java reads it twice): its PASS records are split into SNVs and INDELs. Each
SNV's trinucleotide context and each INDEL's flanking bases (homology and repeats) come from the
reference; SVs are taken once per breakend pair and binned by type and length. The prediction merges the
counts into the model's 12 features (fractions of SNVs by substitution, of INDELs by context, of SVs by
type and length), walks the 500 trees of the forest and writes R's numbers: the probabilities are vote
fractions, formatted as write.table formats doubles (15 significant digits). The forest is
`resources/chord_forest.tsv`, exported once from CHORD.rds by `scripts/export_forest.R` (REPORT.md
section 2) and embedded in the binary.

## Tests

`cargo test --release` runs the Java unit tests of CHORD, ported (`src/tests/`, one module per Java test
class, with the Java test resources in `testdata/`; ChordApplicationTest, which Java ignores as it needs R,
runs here), and tests of the port (`src/tests/port.rs`: each corrected finding in both modes, the QC
rules, the forest read from CHORD.rds against the embedded one, the HG008 prediction;
`src/predict/rfmt.rs`: R's number formatting against 1,742 values written by R, R's sample id conversion,
R's Mersenne-Twister). The comparisons with Java are in [REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
