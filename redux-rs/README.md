# redux-rs

A Rust port of [REDUX](../redux/README.md) 2.1. It takes the same command line and writes the same
files: the duplicate-marked BAM with consensus reads, `<sample>.redux.duplicate_freq.tsv`,
`<sample>.redux.bqr.tsv`, `<sample>.redux.jitter_params.tsv` and `<sample>.redux.ms_table.tsv.gz`.

By default the port corrects the faults found in the Java code while porting it. With `-java_compat`
it reproduces them, and its output is then identical to Java's; that is how the port was verified.
See "Corrections and `-java_compat`" below.

Status and measurements are in [PLAN.md](PLAN.md). What has been compared with Java, and what Java does
that is worth knowing, is in [REPORT.md](REPORT.md).

## Build

```
cargo build --release
```

The default build uses the mimalloc allocator, which is about 13% faster here than the system one;
`cargo build --release --no-default-features` builds without it.

The binary is `target/release/redux-rs`. It links htslib statically (through `rust-htslib`), so samtools
is not needed at run time.

## Run

```
redux-rs \
    -sample SAMPLE_ID \
    -input_bam lane1.bam,lane2.bam \
    -ref_genome ref.fasta -ref_genome_version 38 \
    -unmap_regions unmap_regions.38.tsv \
    -ref_genome_msi_file msi_jitter_sites.38.tsv.gz \
    -form_consensus \
    -output_dir out/ -threads 16
```

Options have the names and defaults of the Java tool. `-bamtool` is accepted and ignored: sorting,
concatenation and indexing are done in the process.

## Corrections and `-java_compat`

The faults found in Java REDUX 2.1 are numbered R1 to R40 in REPORT.md, section 3. For 24 of them the
port takes a corrected path by default; section 3.5 says what each correction does and what it
changes.

- `-java_compat` (or `REDUX_JAVA_COMPAT=1` in the environment) gives Java's behaviour for every
  finding. Use it to compare with Java, or to continue an analysis made with Java REDUX.
- `-java_compat_except R7,R2` (or `REDUX_JAVA_COMPAT_EXCEPT`) gives Java's behaviour for every finding
  but those listed.
- What the corrections change on real data: on Illumina input the BAM stays identical to Java's, apart
  from regions deep enough to hold 100,000 reads at once (R14, R21); on SBX input a few records in ten
  thousand change (R6: a duplicate that shared its bases with the consensus read; R7: forward and
  reverse reads collapsed together). The side files change more: `bqr.tsv` (R1, R2: positions dropped
  or cut off at the end of a 100K block), `jitter_params.tsv` and `ms_table.tsv.gz` (R11:
  microsatellite sites with many rejected reads are now left out) and the labels of
  `duplicate_freq.tsv` (R18). For R7 and R14 the intent of the Java code is not certain; the readings
  taken are in REPORT.md and are open to correction by the REDUX authors.
- Runs that end with an exception in Java complete: `-umi_base_diff_stats` (R24), duplex UMIs without
  the delimiter (R25), Ultima with `-ultima_drop_qual_tags` and the recalibration (R26), the MSI
  prediction for Ultima with the supplied model files (R28). Ultima with `-form_consensus` still stops
  where the reads of a duplicate group disagree: Java has no rule for it (R23).
- The corrected output cannot be checked against Java. It is checked by unit tests of the corrected
  paths, by measuring what each correction changes, and by comparing runs at different thread counts.
- `REDUX_WORKERS=n` in the environment limits the worker threads to `n` without changing how
  `-threads` splits the genome into partitions. It is for comparing with a Java run made with more
  threads than the machine should be given.
- `redux-rs diff_bams a.bam b.bam [threads]` counts the differences between two BAMs of the same reads
  by kind (duplicate flag, cigar, bases, qualities, attributes, records present in one only).

## What "the same" means

All of this section is about the port run with `-java_compat`.

- **BAM:** the same records in the same order with the same fields and attributes, compared after
  decompression. The compressed bytes differ, since the compressor differs.
- **Other files:** byte-identical.
- **Thread count:** Java's result depends on the thread count in three ways; the port follows the first
  and removes the other two:
  - The genome is split into `2 x threads` partitions and duplicates are not matched across a partition
    boundary. The port splits the genome the same way, so it matches Java at the same thread count.
  - With more than one thread, Java writes the unplaced (fully unmapped) reads at the end of the BAM in
    an order that depends on thread timing. The port always writes them in the order Java uses with one
    thread.
  - In regions deep enough to hold 100,000 reads at once, Java's duplicate groups depend on which
    partitions the same thread processed earlier (REPORT.md, section 2). The port always gives the
    result of one thread taking the partitions in order.

  The default mode also removes the third at its root (R21: the cache state starts again with each
  region), and `bqr.tsv` and `duplicate_freq.tsv` no longer depend on the thread count (R16, R18).
  What is left is the first: a duplicate group cut by a partition boundary. On the test inputs that
  was 3 records of 27 million at 128 partitions and none at 16 (REPORT.md, section 2).

## Differences from the Java tool

- **No external tool.** Java calls samtools to sort the interim BAM of unmapped reads, to concatenate the
  per-partition BAMs and to index the result. The port sorts the reads it moves in small per-partition files,
  joins the partition files by copying compressed blocks, and builds the index while writing.
- **Memory.** A few hundred MB instead of about 20 GB.
- **`-compress_level N`** (not a Java option): the deflate level of the output BAM, 0 to 9, default 4.
  Level 4 gives a file within about 0.5% of Java's size (smaller on some inputs, larger on others);
  level 5 is about 2.5% smaller and 7% slower, level 3 about 5% larger and 9% faster, level 1 about 7%
  larger and 21% faster. Compression is about half of the CPU time.
- **`compare_bams a.bam b.bam [threads]`** as the first argument compares two BAMs record by record. It
  treats integer attributes of different widths as equal and accepts a different order of the unplaced
  reads at the end. **`diff_bams`** with the same arguments counts the differences by kind instead of
  stopping at the first.
- **`-java_compat`, `-java_compat_except <ids>`** (not Java options): see above.
- Log lines have the Java layout without the thread name. Debug-only options (`-run_checks`,
  `-log_read_ids`, `-perf_log_time`, `-log_dup_group_size`, `-log_read_cache`, `-read_output`) are
  accepted; see PLAN.md for which have an effect.
