# bam-tools-rs

BamTools metrics in Rust: a port of hmftools BamTools 1.6.2 `com.hartwig.hmftools.bamtools.metrics.BamMetrics`,
the BamTools command oncoanalyser runs on every tumor and normal BAM. It reads the BAM (or CRAM) once and
writes the same files as the Java tool: coverage, flag, fragment-length and per-partition statistics,
and the depth over the driver panel's genes and exons. With `-java_compat` the output is identical to
Java's; by default the faults found in the Java code are corrected (../NOTES.md, BT1-BT20).

## Running

```
bam-tools-rs [com.hartwig.hmftools.bamtools.metrics.BamMetrics | metrics] <options>
```

The command may be left out (the Java jar's main class is BamMetrics). JVM options (`-Xmx...`, `-XX:...`,
`-D...=...`) are ignored; `rust-docker/bin/bamtools` is a wrapper with the Java tool's name for the pipeline.
`bam-tools-rs -version` prints `BamTools version 1.6.2`; `-help` lists the options.

The oncoanalyser fork runs Hartwig's `bam-tools:1.7-beta.3` image, built from the tag
`bam-tools-v1.7-beta.3`, whose BamTools and hmf-common metrics sources are the ones this port follows
(the repository's pom still says 1.6.2).

oncoanalyser's call:

```
bam-tools-rs com.hartwig.hmftools.bamtools.metrics.BamMetrics -sample S -bam_file S.bam -ref_genome ref.fasta \
  -ref_genome_version 38 -driver_gene_panel DriverGenePanel.38.tsv -ensembl_data_dir ensembl_data \
  -threads 2 -output_dir out/   [-regions_file panel.bed]
```

## Options

| Option | Default | Meaning |
|---|---|---|
| `-bam_file` | required | BAM or CRAM, indexed |
| `-ref_genome` | required | reference FASTA (used to decode a CRAM) |
| `-ref_genome_version` | from the BAM header | 37 or 38 (V37, V38, HG37, HG38) |
| `-sample` | "null" | sample id, the prefix of the output files |
| `-output_dir` | the BAM's directory | output directory |
| `-output_id` | | suffix of the partition stats and off-target files |
| `-threads` | 1 | worker threads; 0 or 1 for one |
| `-partition_size` | 1000000 | bases per partition |
| `-map_qual_threshold` | 20 | reads below it count as LOW_MAP_QUAL and add no depth |
| `-base_qual_threshold` | 10 | bases below it count as LOW_BASE_QUAL |
| `-max_coverage` | 250 (10000 with `-regions_file`) | depth cap; bases above it count as MAX_COVERAGE |
| `-exclude_zero_coverage` | | leave zero-depth bases out of the distribution and statistics |
| `-no_duplicates` | | the BAM has its duplicates removed (consensus reads count as reads) |
| `-specific_chr`, `-specific_regions` | | restrict to chromosomes (`chr1;chr2`) or regions (`chr1:100-200;...`) |
| `-regions_file` | | target regions (BED, or TSV/CSV with Chromosome/PositionStart/PositionEnd); statistics over these bases, target region fragment counts |
| `-only_target` | | read only the target regions |
| `-write_off_target` | | off-target fragment pile-ups (overlap counts) |
| `-off_target_frag_overlap_threshold` | 0 | write off-target pile-ups at least this deep |
| `-driver_gene_panel`, `-ensembl_data_dir` | | both needed for the gene and exon coverage |
| `-partition_read_count_check` | 1000000 | debug log every n reads of a partition |
| `-log_read_ids` | | debug log of these reads (';'-separated) |
| `-log_debug`, `-log_level` | | log level |
| `-perf_debug` | | timing summary |
| `-java_compat`, `-java_compat_except BT1,...` | | Java's behaviour for every finding, or for all but those listed; also `BAMTOOLS_JAVA_COMPAT=1`, `BAMTOOLS_JAVA_COMPAT_EXCEPT` |
| `-file_buffer_kb` | 1024 | port only: size of htslib's file reads in KB (0 for htslib's default) |

## Outputs

In `-output_dir`, named `<sample>.bam_metric.<type>.tsv`:

| Type | Content |
|---|---|
| `summary` | territory, read counts, mean/sd/median/MAD depth, filtered base percentages, share of bases at 1x-100x |
| `coverage` | bases per depth (collapsed into ~1% buckets for a cap above 3000) |
| `frag_length` | fragments per length |
| `flag_counts` | samtools-flagstat style counts |
| `partition_stats[.<id>]` | reads per partition, with the processing time (not written with `-only_target`) |
| `gene_coverage`, `exon_coverage` | depth distribution of each panel gene and exon (canonical coding exons with 10-base splice margins); exons with median and mean depth |
| `target_regions` | with `-regions_file`: fragments per target region |
| `off_target_overlaps[.<id>]` | with `-write_off_target`: off-target fragments by overlap count |
| `off_target_enriched[.<id>]` | with `-regions_file` and a threshold: deep off-target pile-ups |

The partition stats and target region rows are in partition order (Java writes them as its threads
finish). `ProcessTime` is a timing and differs between runs.

## Building and testing

```
CARGO_BUILD_JOBS=2 cargo build --release        # target/release/bam-tools-rs
CARGO_BUILD_JOBS=2 cargo test --release         # 23 tests: the Java unit tests ported, the findings,
                                                # and 704 synthetic cases recorded from the Java classes
```

`testdata/java_harness.tsv` comes from `harness/BtHarness.java` (compiled against the BamTools jar);
The comparisons with Java are in [REPORT.md](REPORT.md).
It depends on `../hmf-common-rs` (the Ensembl cache, the driver panel and the metrics file formats live
there). Build with `--no-default-features` for the system allocator instead of mimalloc.

## Speed

On the small tumor/normal pair at 2 threads it takes 12-15 s and 77 MB where Java takes 106-119 s and
2.2 GB; a 100x whole genome takes about 8-9 minutes and 110 MB (../NOTES.md, section 5). Most of the
time is BGZF decompression.
