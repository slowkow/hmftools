# sage-rs

A Rust reimplementation of [SAGE](../sage/README.md) 5.1, the somatic SNV / MNV / indel caller.
It takes the same command-line options and aims to write the same VCF, faster and with less memory.

See [PLAN.md](PLAN.md) for status and measurements and [QUIRKS.md](QUIRKS.md) for Java behaviours that
are reproduced deliberately.

## Supported

Tumor/normal somatic calling on Illumina BAMs with REDUX outputs
(`<sample>.redux.bqr.tsv`, `.redux.jitter_params.tsv`, `.redux.ms_table.tsv.gz`).
Tumor-only runs are implemented but have not been compared with Java.
Not yet: TINC, append mode, CRAM, germline mode, `high_depth_mode`, more than one tumor or reference
sample, SBX, Ultima, the visualiser. Unsupported options stop with an error rather than being ignored.

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
`-log_evidence` writes the per-sample genotype fields of every candidate that passes the hard filters.

## Comparing with Java

```
cp bench/env.example.sh bench/env.local.sh     # edit paths; not committed
bench/run_java.sh chr21 1 -specific_chr chr21
bench/run_rust.sh chr21 1 -specific_chr chr21
diff <(gzip -dc $SAGE_OUT_DIR/java/chr21/chr21.sage.vcf.gz) <(gzip -dc $SAGE_OUT_DIR/rust/chr21/chr21.sage.vcf.gz)
```

Run Java with one thread when comparing: its phase set ids depend on thread scheduling.
