# cobalt-rs

A Rust port of hmftools COBALT 3.1 (`../cobalt`): the read depth of a tumor and/or reference BAM or CRAM
in 1,000-base windows, GC-normalised into ratios, and their segmentation (the `.pcf` files PURPLE reads).
It takes the same options and writes the same files as the Java `CobaltApplication`; with `-java_compat`
the files are byte-identical to Java's (the gzip bytes too). Faults found in the Java code are corrected by
default (findings CO1-CO6 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/cobalt-rs
cobalt-rs -tumor TUMOR -tumor_bam tumor.bam -reference NORMAL -reference_bam normal.bam \
    -gc_profile GC_profile.1000bp.38.cnp -ref_genome_version 38 -threads 8 -output_dir cobalt/
```

Modes: tumor/normal; tumor-only (no reference; whole genome needs `-tumor_only_diploid_bed`); germline-only
(no tumor); targeted (`-target_region_norm_file`, one or more panels by the `PanelLabel` column).

JVM options (`-Xmx16g`, `-XX:...`) are accepted and ignored, so the bioconda wrapper's command line works;
`cobalt-rs -version` prints `Cobalt version 3.1`.

## Options

| Option | Default | |
|---|---|---|
| `-tumor`, `-tumor_bam` | | tumor id and BAM/CRAM |
| `-reference`, `-reference_bam` | | reference id and BAM/CRAM |
| `-gc_profile` | required | GC profile (uncompressed) |
| `-ref_genome_version` | required | 37 or 38 (V37, HG37, ...) |
| `-ref_genome` | | FASTA, for CRAM |
| `-tumor_only_diploid_bed` | | diploid regions (tumor-only mode) |
| `-target_region_norm_file` | | targeted panel normalisation file |
| `-min_quality` | 10 | minimum mapping quality |
| `-include_duplicates` | | count duplicates (and drop consensus reads, attribute CR) |
| `-gc_ratio_min`, `-gc_ratio_max` | 0.24, 0.68 | GC limits (no effect in Java, CO1) |
| `-pcf_gamma` | 100 | segmentation penalty factor |
| `-skip_pcf_calc` | | no .pcf files |
| `-specific_chr`, `-specific_regions` | | restrict to these chromosomes (whole chromosomes, as Java, CO6) |
| `-threads` | 1 | worker threads (0 is taken as 1; Java fails, CO5) |
| `-output_dir` | | created if missing |
| `-bam_validation` | SILENT | accepted, not applied |
| `-log_debug`, `-log_level` | | debug logging, with `PerfStats` timings |
| `-java_compat`, `-java_compat_except CO1,...` | | Java behaviour for every finding, or all but those listed (also `COBALT_JAVA_COMPAT=1`, `COBALT_JAVA_COMPAT_EXCEPT`) |

## Outputs

In `-output_dir`:
- `<tumor or reference>.cobalt.ratio.tsv.gz`: per window chromosome, position, read depths, GC ratios,
  the reference's diploid-adjusted ratio and the GC contents (-1 for masked or absent values)
- `<sample>.cobalt.ratio.pcf`, tumor and reference: the segments (`Chromosome Start End MeanRatio`)
- `<tumor>.cobalt.gc.median.tsv`, `<reference>.cobalt.gc.median.tsv`: mean and median depth, smoothed
  median per GC bucket
- `<reference>.cobalt.ratio.median.tsv`: median ratio and window count per chromosome
- `cobalt.version`

## How it works, and speed

Each chromosome is read in 10 Mb partitions by `-threads` workers (tumor and reference share them), each
with its own htslib reader, counting the aligned bases and G/C bases of passing reads per window into its
own slice of the window arrays; the result does not depend on threads or partitioning. Then the GC buckets
and their smoothed medians, the GC, mean/median, low-coverage, diploid and unity normalisations, and the
segmentation of hmf-common-rs. Reading runs at about the speed of `samtools view -c`; for CRAM only the
needed fields are decoded. Memory is a few hundred MB (Java: GB).

## Tests

```
cargo test --release
```

runs the ported Java unit tests (`src/java_tests.rs`), own tests, and the Java e2e scenarios
(`src/e2e_tests.rs`) on the Java tests' synthetic BAMs, found in `../cobalt/src/test/resources` or where
`COBALT_TEST_RESOURCES` points (skipped otherwise). `COBALT_E2E_KEEP=<dir>` keeps each scenario's inputs
and arguments for a comparison with Java. The comparisons with Java are in [REPORT.md](REPORT.md).
