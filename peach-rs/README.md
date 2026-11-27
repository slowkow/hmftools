# peach-rs

A Rust port of hmftools PEACH 2.0 (`../peach`): the pharmacogenomic haplotypes of a panel's genes (DPYD and
UGT1A1 in the resource panel) called from a sample's germline variants, with each haplotype's function and
the drugs and prescription information linked to the gene. It takes the same options and writes the same
five files as the Java `PeachApplication`; with `-java_compat` they are byte-identical to Java's. Faults
found in the Java code are corrected by default (findings PE1-PE2 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/peach-rs
peach-rs -sample_name NORMAL -vcf_file TUMOR.purple.germline.vcf.gz -haplotypes_file haplotypes.38.tsv \
    -function_file haplotype_functions.38.tsv -drugs_file peach_drugs.38.tsv -log_level INFO -output_dir peach/
```

This is the command line oncoanalyser runs (PEACH). JVM options are accepted and ignored, and
`peach-rs -version` prints `Peach version 2.0`.

## Options

`-vcf_file` (bgzipped or plain), `-sample_name`, `-haplotypes_file`, `-drugs_file` and `-function_file`
(both optional), `-output_dir` (created if missing), `-log_level`, `-log_debug`; `-java_compat` and
`-java_compat_except PE1,...` (or `PEACH_JAVA_COMPAT=1` and `PEACH_JAVA_COMPAT_EXCEPT`).

## Outputs

In `-output_dir`, each sorted as Java sorts it:
- `<sample>.peach.events.tsv`: each event of the VCF at a panel position with its count (0, 1, 2 or
  UNKNOWN for a no-call);
- `<sample>.peach.gene.events.tsv`: the events per gene;
- `<sample>.peach.haplotypes.all.tsv`: every combination of haplotypes that explains a gene's events;
- `<sample>.peach.haplotypes.best.tsv`: the best combination's haplotypes with function, linked drugs and
  prescription URLs, or "Unresolved Haplotype";
- `<sample>.peach.qc.tsv`: the QC status per gene.

## How it works

An event is a VCF record (not filtered) whose reference bases cover a position of a panel haplotype's
event. Per gene, the haplotype combinations that explain the event counts are found by recursive descent
over the gene's haplotypes (as HaplotypeCaller), the default haplotype filling up to two copies; the best are
those nearest to two haplotypes, then with the fewest haplotypes other than the wild type. Every output is
sorted, so nothing depends on Java's HashMap order.

## Tests

`cargo test --release` runs the Java unit tests (PeachAlgoTest with TestHaplotypePanelFactory,
PeachConfigTest, the loader tests, DrugInfoStoreTest, HaplotypeCombinationTest, HaplotypePanelTest and the
five output tests; `src/tests/java_algo.rs` and `src/tests/java_io.rs`, with the Java tests' resource files in
`testdata/`) and tests of the port (`src/tests/port.rs`: event ids, chromosome names, PE1 in both modes).
