# sage-rs

A Rust reimplementation of [SAGE](../sage/README.md) 5.1, the somatic SNV / MNV / indel caller.
It takes the same command-line options and aims to write the same VCF, faster and with less memory.

See [PLAN.md](PLAN.md) for status and measurements and [QUIRKS.md](QUIRKS.md) for Java behaviours that
are reproduced deliberately.

## Supported

Tumor/normal somatic calling on Illumina, SBX and Ultima (`-sequencing_type SBX` / `ULTIMA`) BAMs or
CRAMs (versions 3.0 and 3.1, decoded against `-ref_genome`) with REDUX outputs
(`<sample>.redux.bqr.tsv`, `.redux.jitter_params.tsv`, `.redux.ms_table.tsv.gz`).
Also: tumor-only runs, several tumor and reference samples, `-panel_only`, `-high_depth_mode`,
`-germline`, samples REDUX predicts as MSI, split (RNA) reads, TINC (`-run_tinc` with `-pon_file` and
`-gnomad_freq_dir`, or standalone as `sage-rs tinc -input_vcf ...`; with several tumors it uses the
first, as Java does), append mode
(`sage-rs append -input_vcf ... -reference ... -reference_bam ... -output_vcf ...`),
`-write_frag_lengths` (`<output>.frag_lengths.tsv.gz`) and `-specific_var_file` (a delimited file with
`Chromosome`, `Position`, `Ref` and `Alt` columns; ignored, as in Java, when `-specific_chr` or
`-specific_regions` is given).
Not yet: the visualiser.

CRAM input gives slightly different calls from BAM input of the same reads, in Java and so here too:
see "CRAM" in [QUIRKS.md](QUIRKS.md).

## Build and run

```
cargo build --release
target/release/sage-rs \
    -tumor TUMOR -tumor_bam tumor.redux.bam \
    -reference NORMAL -reference_bam normal.redux.bam \
    -jitter_bqr_dir /dir/with/redux/files \
    -ref_genome ref.fasta -ref_genome_version 38 \
    -hotspots KnownHotspots.somatic.38.vcf.gz \
    -driver_gene_panel DriverGenePanel.38.tsv \
    -high_confidence_bed highconf.bed.gz \
    -ensembl_data_dir ensembl_data/ \
    -include_mt -threads 16 \
    -output_vcf out.sage.vcf.gz
```

Extra options: `-log_candidates` writes `<output>.candidates.tsv` (same columns as Java);
`-log_evidence` writes the per-sample genotype fields of every candidate that passes the hard filters;
`-read_cache_mb N` sets how much tumor read data each thread keeps between the candidate and evidence
stages (default 32; 0 reads the file twice, as Java does, for about 0.5 GB less memory at 16 threads
and 15% more time); `-perf_summary` reports where the wall-clock time went.

## Comparing with Java

```
cp bench/env.example.sh bench/env.local.sh     # edit paths; not committed
bench/run_java.sh chr21 1 -specific_chr chr21
bench/run_rust.sh chr21 1 -specific_chr chr21
diff <(gzip -dc $SAGE_OUT_DIR/java/chr21/chr21.sage.vcf.gz) <(gzip -dc $SAGE_OUT_DIR/rust/chr21/chr21.sage.vcf.gz)
```

Run Java with one thread when comparing: its phase set ids depend on thread scheduling.
`bench/compare_vcf.py` compares two VCFs field by field and reports phase-set-id-only differences
separately, for runs where Java used more than one thread. `SAGE_TUMOR_ONLY=1` makes both run scripts
leave out the reference sample.

Java's partition size depends on the thread count for contigs shorter than `threads` x 100 kb (chrM),
and so do its calls there; compare those with one thread on both sides.

CRAM: `bench/make_cram_set.sh <out dir> [region] [version=3.0]` converts the sample set's BAMs and
writes an `env.sh` for them; then `SAGE_ENV=<out dir>/env.sh SAGE_BAM_EXT=cram bench/run_java.sh ...`
and the same for `run_rust.sh`. The htsjdk in Java SAGE 5.1 reads CRAM 3.0 only.

TINC: `bench/run_tinc_pair.sh <name> <sage.vcf.gz>` runs the standalone Java and Rust analysers on the
same VCF and compares the results; `bench/simulate_tinc.py` adds simulated contamination to a VCF so
that the recovery path is exercised.

## Tests

`cargo test --release` runs the unit tests, including the Java SAGE unit tests ported to `src/tests`
(one module per Java test class). Tests that are not from Java cover the port's own alignment reading
(`alignment_reading`: region queries, the tag parser, the read cache against the file, CRAM against
BAM) and each fast path against the general code it stands in for.

`SAGE_ENV=/path/to/other-env.sh` points the run scripts at another sample set, for example the public
HG008 SBX slices.
