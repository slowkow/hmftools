# sage-rs

A Rust reimplementation of [SAGE](../sage/README.md) 5.1, the somatic SNV / MNV / indel caller.
It takes the same command-line options and aims to write the same VCF, faster and with less memory.

See [PLAN.md](PLAN.md) for status and measurements and [QUIRKS.md](QUIRKS.md) for Java behaviours that
are reproduced deliberately.

## Supported

Tumor/normal somatic calling on Illumina, SBX and Ultima (`-sequencing_type SBX` / `ULTIMA`) BAMs with REDUX outputs
(`<sample>.redux.bqr.tsv`, `.redux.jitter_params.tsv`, `.redux.ms_table.tsv.gz`).
Also: tumor-only runs, several tumor and reference samples, `-panel_only`, `-high_depth_mode`,
`-germline`, samples REDUX predicts as MSI, split (RNA) reads, TINC (`-run_tinc` with `-pon_file` and
`-gnomad_freq_dir`, or standalone as `sage-rs tinc -input_vcf ...`) and append mode
(`sage-rs append -input_vcf ... -reference ... -reference_bam ... -output_vcf ...`).
Not yet: CRAM, the visualiser, `-write_frag_lengths`, `-specific_var_file`.

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
`bench/compare_vcf.py` compares two VCFs field by field and reports phase-set-id-only differences
separately, for runs where Java used more than one thread. `SAGE_TUMOR_ONLY=1` makes both run scripts
leave out the reference sample.

TINC: `bench/run_tinc_pair.sh <name> <sage.vcf.gz>` runs the standalone Java and Rust analysers on the
same VCF and compares the results; `bench/simulate_tinc.py` adds simulated contamination to a VCF so
that the recovery path is exercised.

## Tests

`cargo test --release` runs the unit tests, including the Java SAGE unit tests ported to `src/tests`
(one module per Java test class).

`SAGE_ENV=/path/to/other-env.sh` points the run scripts at another sample set, for example the public
HG008 SBX slices.
