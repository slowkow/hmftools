# redux-rs

A Rust port of [REDUX](../redux/README.md) 2.1. It takes the same command line and writes the same
files: the duplicate-marked BAM with consensus reads, `<sample>.redux.duplicate_freq.tsv`,
`<sample>.redux.bqr.tsv`, `<sample>.redux.jitter_params.tsv` and `<sample>.redux.ms_table.tsv.gz`.

Status and measurements are in [PLAN.md](PLAN.md). What has been compared with Java, and what Java does
that is worth knowing, is in [REPORT.md](REPORT.md).

## Build

```
cargo build --release
```

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

## What "the same" means

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

## Differences from the Java tool

- **No external tool.** Java calls samtools to sort the interim BAM of unmapped reads, to concatenate the
  per-partition BAMs and to index the result. The port sorts the reads it moves in small per-partition files,
  joins the partition files by copying compressed blocks, and builds the index while writing.
- **Memory.** A few hundred MB instead of about 20 GB.
- **`-compress_level N`** (not a Java option): the deflate level of the output BAM, 1 to 12, default 4.
  Level 4 gives a file within about 0.5% of Java's size (smaller on some inputs, larger on others);
  level 5 is about 2.5% smaller and 7% slower.
- **`compare_bams a.bam b.bam [threads]`** as the first argument compares two BAMs record by record. It
  treats integer attributes of different widths as equal and accepts a different order of the unplaced
  reads at the end.
- Log lines have the Java layout without the thread name. Debug-only options (`-run_checks`,
  `-log_read_ids`, `-perf_log_time`, `-log_dup_group_size`, `-log_read_cache`, `-read_output`) are
  accepted; see PLAN.md for which have an effect.
