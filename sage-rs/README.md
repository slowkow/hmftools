# sage-rs

A Rust reimplementation of [SAGE](../sage/README.md) 5.1, the somatic SNV / MNV / indel caller.
It takes the same command-line options, and is faster and uses less memory.

The Java reference is SAGE 5.1 as of hmftools commit `a224a2cd8e`. The Java source beside this crate is
newer: upstream commit `de50386378` (2026-11) changed SAGE without changing its version number (the
realigned-support test for indels near a soft-clipped end, a base-quality penalty scaled between 6
and 12, SAGE's own `NM` count, and append's version check). Those changes are not in the port yet, so
`-java_compat` reproduces the earlier commit, not the current Java source.

By default the port corrects the faults found in the Java code. With `-java_compat` it behaves as Java
does and writes the same VCF, byte for byte; that is how the port was verified. See "Corrections and
`-java_compat`" below.

See [PLAN.md](PLAN.md) for status and measurements, [REPORT.md](REPORT.md) for what was compared with
Java and for the numbered list of Java findings (section 5), and [QUIRKS.md](QUIRKS.md) for the same
findings by subject.

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
`-specific_regions` is given), and the visualiser (below).

CRAM and BAM input of the same reads give the same calls. In Java, and here with `-java_compat`, they
differ slightly (finding S7 in [REPORT.md](REPORT.md)).

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
and 15% more time); `-file_buffer_kb N` sets the buffer each alignment file handle reads through
(default 1024; 0 leaves htslib's default, which is about half as fast on a whole genome);
`-perf_summary` reports where the wall-clock time went.

**Base-quality bands** (options of the port, not of Java SAGE). Java fixes two quality bands for
Illumina data: a base at Q26 or below is low quality (a read with such a base at the variant cannot
create a candidate by itself), and only a base at Q30 or above is high quality (`minTumorSupport` needs
3 strong-support reads of high quality, 2 at a hotspot; high-depth mode drops the others as evidence).
On data whose base qualities were rewritten by recalibration and stop below Q30, no SNV passes.
`-low_base_qual_threshold N` and `-high_base_qual_threshold N` move the bands (defaults 26 and 30; the
low band must end below the high one). Not given, the output is unchanged, with or without
`-java_compat`; given, the VCF header records the bands (`##sageBaseQualThresholds=low:N,high:N`) and
the log says so. The `minAvgRawBaseQual` filter has its own Java options (`-min_avg_base_qual`, 25, and
`-min_avg_base_qual_hotspot`, 18). Hartwig tuned its filters on data with the usual qualities, so calls
made with moved bands need checking, for example against another sample of the same tumor.

`cargo build --release --features mimalloc` links the mimalloc allocator: the same output for about 16%
less CPU time, at several times the memory (11.6 GB against 1.7 GB on a whole genome). It only pays
when storage is fast enough to keep every thread busy.

## Corrections and `-java_compat`

The faults found in Java SAGE are numbered S1 to S48 in [REPORT.md](REPORT.md), section 5. For 29 of
them the port takes a corrected path by default; section 5.5 says what each correction does and what
it changes on the test data.

- `-java_compat` (or `SAGE_JAVA_COMPAT=1` in the environment) gives Java's behaviour for every finding.
  Use it to compare with Java, or to continue an analysis made with Java SAGE. It works for the caller,
  `sage-rs append` and `sage-rs tinc`.
- `-java_compat_except S4,S7` (or `SAGE_JAVA_COMPAT_EXCEPT`) gives Java's behaviour for every finding
  but those listed.
- The corrections that change the most records are the three in local phasing (S1 to S3: which variants
  share an `LPS` id), the CRAM `NM` tag (S7), and in `append` the formatting of the existing records
  (S19). For S5 (a pair of mates whose merge fails), S10 (the early exit of the TINC fit) and S45 (the
  `AED` row of the visualiser) the intent of the Java code is not certain; the readings taken are
  described in REPORT.md and are open to correction by the SAGE authors.
- Twelve findings concern Java output that is not reproducible or a Java failure (thread-order `LPS`
  ids, the `NEAR_HOTSPOT` flag, the random choice of reads in the visualiser, CRAM 3.1, ...). There
  the port does one fixed thing in both modes, as before.
- The corrected output cannot be checked against Java. It is checked by unit tests of the corrected
  paths (`src/tests/corrections.rs` and the tests named in REPORT.md 5.5), by being the same at any
  thread count, and by measuring what each correction changes.
- The header of a corrected VCF differs from Java's in two declarations: `RC_JIT` is declared with the
  two values that are written (S40), and the `fragmentLength` filter, which nothing applies, is not
  declared (S41).

## Visualiser

As Java SAGE, `sage-rs` can write one gzipped HTML file per variant showing the reads behind it
(`<sample><chromosome>_<position>_<ref>_<alt>.html.gz`):

| Option | Meaning |
|---|---|
| `-vis_variants chr:pos:ref:alt;...` | Write files for these variants. Each also adds the 200 bases either side of it to the regions that are called, so the run is limited to those regions. Give the variants of a chromosome in position order: phase set ids follow the order given, and Java can fail on variants that are out of order. |
| `-vis_pass_only` | Write files for passing variants only (for every passing variant when no variants are named). |
| `-vis_max_support_reads N` | Reads shown per match type: 0 (default) for the built-in limits (40 full, 20 ref and realigned, 10 core and partial core, 5 others), `-1` for up to 1000, otherwise `N`. |
| `-vis_output_dir DIR` | Where the files go; by default `vis/` beside the output VCF. |
| `-vis_purple_vcf FILE` | Write files for the variants this PURPLE VCF flags `REPORTED`, with an amino acid track drawn from their `IMPACT` annotation and `-ensembl_data_dir`; the file names then carry the gene and the coding effect. |

With `-java_compat` the files are the same as Java's apart from which reads are shown when there are
more than the limit and the order of reads that tie, both of which Java leaves to a random number
generator, and text positions, which in Java depend on the fonts of the machine: see
[QUIRKS.md](QUIRKS.md). By default six faults of the Java visualiser are corrected (S38, S39, S44 to
S47): among them, a reportable variant whose amino acid track cannot be drawn still gets its file,
without the track. Append mode has no visualiser, as in Java.

## Comparing with Java

Everything in this section needs the port in Java mode: `export SAGE_JAVA_COMPAT=1` before
`bench/run_rust.sh`, or add `-java_compat` to its arguments.

```
cp bench/env.example.sh bench/env.local.sh     # edit paths; not committed
export SAGE_JAVA_COMPAT=1
bench/run_java.sh chr21 1 -specific_chr chr21
bench/run_rust.sh chr21 1 -specific_chr chr21
diff <(gzip -dc $SAGE_OUT_DIR/java/chr21/chr21.sage.vcf.gz) <(gzip -dc $SAGE_OUT_DIR/rust/chr21/chr21.sage.vcf.gz)
```

Run Java with one thread when comparing: its phase set ids depend on thread scheduling.
`bench/compare_vcf.py` compares two VCFs field by field and reports phase-set-id-only differences
separately, for runs where Java used more than one thread. `SAGE_TUMOR_ONLY=1` makes both run scripts
leave out the reference sample.

Java's partition size depends on the thread count for contigs shorter than `threads` x 100 kb (chrM),
and so do its calls there (S24); compare those with one thread on both sides. Without `-java_compat`
the port cuts such a contig the same way at every thread count.

CRAM: `bench/make_cram_set.sh <out dir> [region] [version=3.0]` converts the sample set's BAMs and
writes an `env.sh` for them; then `SAGE_ENV=<out dir>/env.sh SAGE_BAM_EXT=cram bench/run_java.sh ...`
and the same for `run_rust.sh`. The htsjdk in Java SAGE 5.1 reads CRAM 3.0 only.

TINC: `bench/run_tinc_pair.sh <name> <sage.vcf.gz>` runs the standalone Java and Rust analysers on the
same VCF and compares the results; `bench/simulate_tinc.py` adds simulated contamination to a VCF so
that the recovery path is exercised.

## Tests

`cargo test --release` runs the unit tests, including the Java SAGE unit tests ported to `src/tests`
(one module per Java test class). They run in the default, corrected mode; the few whose expected
values hold a Java fault run inside `compat::with_java(true, ...)`, and `src/tests/corrections.rs`
tests the corrected behaviour beside the Java one. Tests that are not from Java cover the port's own alignment reading
(`alignment_reading`: region queries, the tag parser, the read cache against the file, CRAM against
BAM), each fast path against the general code it stands in for, and the visualiser (`visualiser`:
Java has no tests for it, so the expected values were generated by calling the Java classes from
`src/tests/data/VisFixtures.java`).

`SAGE_ENV=/path/to/other-env.sh` points the run scripts at another sample set, for example the public
HG008 SBX slices.
