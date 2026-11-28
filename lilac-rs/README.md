# lilac-rs

A Rust port of hmftools LILAC 2.1 (`../lilac`): HLA typing of a sample from its BAM or CRAM. For each gene
group it collects the fragments of the HLA genes' coding regions, phases their amino acids, narrows the
alleles of the IPD-IMGT/HLA reference down to candidates, scores every candidate combination ("complex")
by its fragment coverage, allele frequencies and complexity, and reports the winning alleles with their
reference, tumor and RNA coverage, the tumor copy number and somatic variants of each allele, and QC.
It takes the same options and writes the same files as the Java `LilacApplication`; with `-java_compat` the
files are byte-identical to Java's (Java run with `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, see
LI1). Faults found in the Java code are corrected by default (findings LI1-LI25 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/lilac-rs
lilac-rs -sample TUMOR -reference_bam normal.bam -tumor_bam tumor.bam -purple_dir purple/ \
    -ref_genome GRCh38.fasta -ref_genome_version 38 -resource_dir lilac_resources/ \
    -freq_score_penalty 0.0009 -genes ALL -sequencing_type ILLUMINA -threads 2 -output_dir lilac/
```

Modes: tumor/normal (typing from the reference BAM, coverage of the winning alleles in the tumor BAM);
tumor-only (`-tumor_bam` only: typing from the tumor); germline-only (`-reference_bam` only). With
`-purple_dir` (or `-gene_copy_number` and `-somatic_vcf`) the tumor copy number of each allele and the
somatic variants in the HLA coding regions are added; `-rna_bam` adds RNA coverage. Gene groups (`-genes`,
comma separated): `MHC_CLASS_1` (HLA-A, -B, -C, with HLA-Y and HLA-H as pseudogenes), `HLA_DQB1`, `HLA_DPA1`,
`HLA_DPB1`, `HLA_DQA1`, `ALL` (those five), and `HLA_DRB1`, `HLA_DRB` (not in `ALL`). oncoanalyser runs
`-genes ALL -freq_score_penalty 0.0009`, or for panels `-genes MHC_CLASS_1 -freq_score_penalty 0.0018
-targeted_panel`.

JVM options (`-Xmx16g`, `-XX:...`) are accepted and ignored, so the bioconda wrapper's command line works;
`lilac-rs -version` prints `Lilac version 2.1`.

## Options

| Option | Default | |
|---|---|---|
| `-sample` | required | sample id (output file prefix) |
| `-reference_bam`, `-tumor_bam`, `-rna_bam` | | BAM or CRAM files |
| `-sample_data_dir` | | instead of the BAM options: `<dir>/<sample>.hla.bam`, the reference `<sample up to its last T>R.hla.bam`, `<sample>.rna.hla.bam`, PURPLE files in the same folder |
| `-resource_dir` | required | `hla_ref_nucleotide_sequences.csv`, `hla_ref_aminoacid_sequences.csv`, `lilac_allele_frequencies.csv` |
| `-ref_genome`, `-ref_genome_version` | required | FASTA (the HLA region's bases, and CRAM), 37 or 38 |
| `-genes` | ALL | gene groups, as above |
| `-sequencing_type` | ILLUMINA | ILLUMINA, ULTIMA (UQ tag low-qual bases), SBX (medium quals); not ILLUMINA doubles two QC warning thresholds |
| `-targeted_panel` | | doubles the unmatched-haplotype and indel QC warning thresholds |
| `-purple_dir`, `-gene_copy_number`, `-somatic_vcf` | | PURPLE outputs (`<sample>.purple.cnv.gene.tsv`, `<sample>.purple.somatic.vcf.gz`) |
| `-freq_score_penalty` | 0.0009 | weight of the cohort allele frequencies in the score |
| `-min_evidence_factor`, `-min_high_qual_evidence_factor`, `-min_evidence_support`, `-min_depth_filter` | 0.006, 0.003, 1, 10 | support needed for a base or amino acid |
| `-min_fragments_per_allele`, `-min_fragments_to_remove_single` | 7, 40 | phasing |
| `-max_ref_fragments` | 1000 | reference fragments used to rank complexes (0: all) |
| `-top_score_threshold` | 0.005 | complexes kept within this fraction of the top coverage |
| `-hla_y_threshold`, `-fatal_total_low_coverage_positions` | 0.003, 300 | |
| `-actual_alleles`, `-restricted_alleles` | | `;`-separated allele lists (forced solution, restricted reference) |
| `-max_elim_candidates` | 0 | |
| `-write_types` | SUMMARY | also FRAGMENTS, REF_COUNTS, or ALL |
| `-threads` | 1 | 0 is taken as 1 (Java fails, LI7) |
| `-output_dir` | the reference BAM's folder | created if missing |
| `-bam_validation`, `-stack_sample_rate`, `-log_perf`, `-debug_phasing`, `-run_validation` | | accepted (the first three have no effect) |
| `-log_debug`, `-log_level` | | debug logging, with the `LIn:` traces of the findings reached |
| `-java_compat`, `-java_compat_except LI2,...` | | Java behaviour for every finding, or all but those listed (also `LILAC_JAVA_COMPAT=1`, `LILAC_JAVA_COMPAT_EXCEPT`) |

## Outputs

In `-output_dir`, one block of rows per gene group:
- `<sample>.lilac.tsv`: per allele of the solution, reference/tumor/RNA total, unique, shared and wildcard
  fragments, tumor copy number, somatic variant counts by effect
- `<sample>.lilac.qc.tsv`: per gene group the QC status and metrics
- `<sample>.lilac.candidates.coverage.tsv`: the ranked complexes with score, penalties and coverage (Java's
  header line is that of the fragments file, LI2; corrected by default)
- with `-write_types FRAGMENTS`: `<sample>.lilac.fragments.tsv`, `<sample>.lilac.candidates.fragments.tsv`
- with `-write_types REF_COUNTS`: `<sample>.lilac.candidates.aminoacids.txt`, `.candidates.nucleotides.txt`
  (the last gene group's, LI17), `<sample>.lilac.HLA-<gene>.aminoacids.txt` / `.nucleotides.txt`,
  `<sample>.lilac.HLA_Y.coverage.tsv` when HLA-Y is found; in debug mode `HLA_Y.fragments.tsv`

## How it works, and speed

The two reference sequence files are read once (Java reads them again for every gene group) and parsed
per group into token sequences (one `u16` per locus; inserts interned). Each gene's coding region is read
with an indexed htslib query; reads become fragments with the read handling of Java (soft clips, low-qual
trimming, indels matched against the reference's indel alleles by substring search). Phasing,
candidate selection and the fragment-allele mapping follow the Java classes one to one, with Java's HashMap
and HashSet iteration orders reproduced where they can reach the output (`hmf_common::java_hash`, the JDK 19+
and Guava 33 table sizes), TimSort for every sort with a comparator, and fdlibm's `log10`.

The complexes (up to about a million: the product of the per-gene candidate pairs) are scored with exact
integer coverage totals (each shared fraction 1/k is a whole number of 1/27720), summed once per gene
complex or gene pair for the fragments only those genes' alleles support; only the complexes Java keeps get
their per-allele coverage computed in floating point as Java does, and the few whose totals are exactly
half way are recomputed so. Java's culling of low-coverage complexes (per thread, LI6), its
non-transitive ranking comparator and its tie-breaks are reproduced. Gene groups run in parallel with
`-threads 2` or more and are written in order.

A whole-genome tumor/normal pair takes about 0.7 s and 0.3-0.4 GB with 2 threads (Java: 30 s, 5 GB with
-Xmx8g); 982,800 complexes are scored in 0.15 s (Java 5.6 s).

## Tests

```
cargo test --release
```

runs the ported Java unit tests (`src/java_tests.rs`: allele, sequence, fragment, read, coverage, complex,
ranking, evidence, QC and the three application tests on the Java tests' sequences, copied to `testdata/`)
and the port's own tests (Java hash orders, the integer complex totals against floating point on random
fragment sets, the substring search against a port of hmf-common's SuffixTree, fdlibm log10 against the
JVM's values where `/tmp/lilac_hash/log10_java.txt` exists). The comparisons with Java are in
[REPORT.md](REPORT.md); the scripts in `../scripts` (`cases.sh java|rust|list`, `compare.sh`,
`compare_all.sh`, `trace_findings.sh`, `timing_table.sh`) rerun them.

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
