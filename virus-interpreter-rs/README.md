# virus-interpreter-rs

A Rust port of hmftools Virus Interpreter 1.7.1 (`../virus-interpreter`): VIRUSBreakend's summary of the
viruses found in a tumor, annotated with their names, interpretation (EBV, HBV, HHV-8, HPV, MCV) and the
expected clonal coverage, and classified as reported, blacklisted and of high or low driver likelihood. It
takes the same options and writes the same file as the Java `VirusInterpreterApplication`; with
`-java_compat` it is byte-identical to Java's. Faults found in the Java code are corrected by default
(findings VI1-VI3 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/virus-interpreter-rs
virus-interpreter-rs -sample TUMOR -purple_dir purple -tumor_metrics_dir bamtools_TUMOR \
    -virus_breakend_tsv TUMOR.virusbreakend.vcf.summary.tsv -taxonomy_db_tsv taxonomy_db.tsv \
    -virus_reporting_db_tsv virus_reporting_db.tsv -virus_blacklisting_db_tsv virus_blacklisting_db.tsv \
    -log_level INFO -output_dir virusinterpreter/
```

This is the command line oncoanalyser runs (VIRUSINTERPRETER). It names the BamTools folder
`-tumor_metrics_dir`, as Hartwig's virus-interpreter 1.7.2-beta.1 image does; the 1.7.1 jar of this
repository calls it `-bam_metrics_tumor_dir`. The port takes either. JVM options are accepted and ignored,
and `virus-interpreter-rs -version` prints `VirusInterpreter version 1.7.1`.

## Options

`-sample`, `-purple_dir` (`<sample>.purple.qc` and `<sample>.purple.purity.tsv`), `-bam_metrics_tumor_dir` or
`-tumor_metrics_dir` (`<sample>.bam_metric.summary.tsv`), `-virus_breakend_tsv`, `-taxonomy_db_tsv`,
`-virus_reporting_db_tsv`, `-virus_blacklisting_db_tsv`, `-output_dir`, `-log_level`, `-log_debug`;
`-java_compat` and `-java_compat_except VI1,...` (or `VIRUS_INTERPRETER_JAVA_COMPAT=1` and
`VIRUS_INTERPRETER_JAVA_COMPAT_EXCEPT`).

## Output

`<sample>.virus.annotated.tsv`: taxid, name, qcStatus, integrations, interpretation, percentageCovered,
meanCoverage, expectedClonalCoverage (`null` when PURPLE's QC fails for no tumor or contamination),
reported, blacklisted, driverLikelihood; doubles as Java's `Double.toString` writes them.

## How it works

The expected clonal coverage is the tumor's mean coverage times purity over ploidy. A virus of the reporting
database is reported when it passes VIRUSBreakend's QC and, with acceptable PURPLE quality, its covered
percentage and mean depth exceed the database's minimum (integrated or not) and the expected clonal
coverage; without a minimum it is reported regardless, and with unacceptable PURPLE quality only then.
Java holds the 2.6-million-line taxonomy database in memory (about 0.9 GB at peak); the port reads the
summary first and keeps only the names it looks up, while checking every line as Java parses it.

## Tests

`cargo test --release` runs the Java unit tests (VirusInterpreterAlgoTest, VirusReportingDbModelTest,
VirusBlacklistingDbFileTest, VirusReportingDbFileTest, CoveragesAnalyzerTest, TaxonomyDbTest,
TaxonomyDbFileTest; `src/tests/java_tests.rs`, with their resource files in `testdata/`) and tests of the
port (`src/tests/port.rs`).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
