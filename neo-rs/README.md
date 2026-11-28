# neo-rs

A port of hmftools NEO 1.4 (Java) to Rust: the neo-epitopes of a tumor's somatic point mutations and fusions
(NeoEpitopeFinder), and their peptides scored for HLA binding, expression and immunogenic similarity against the
sample's alleles (NeoScorer). These are oncoanalyser's NEO_FINDER and NEO_SCORER steps. The port takes NEO's command
lines and writes the same files; with `-java_compat` they equal Java's byte for byte. Without it, the Java faults
listed in [REPORT.md](REPORT.md) (NE1 to NE14) are corrected.

## Usage

As the bioconda wrapper `neo` takes it: the first argument may be the main class.

```
neo-rs [com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder] -sample <id> -linx_dir <dir> -somatic_vcf <purple vcf> \
    -ref_genome <fasta> -ref_genome_version 38 -ensembl_data_dir <dir> -output_dir <dir>
neo-rs com.hartwig.hmftools.neo.score.NeoScorer -sample <id> -purple_dir <dir> -lilac_dir <dir> -neo_dir <dir> \
    -ensembl_data_dir <dir> -score_file_dir <neo binding dir> -cancer_tpm_medians_file <csv.gz> -output_dir <dir>
```

With no class, or NeoEpitopeFinder (the jar's main class), the finder runs. JVM options (`-Xmx...`, `-Xms...`,
`-XX:...`, `-D...=...`) are accepted and ignored. `neo-rs -version` prints `Neo version 1.4`.

Port options, for either class: `-java_compat` (or `NEO_JAVA_COMPAT=1`) reproduces every Java fault;
`-java_compat_except NE6,NE13` (or `NEO_JAVA_COMPAT_EXCEPT`) all but those listed.

### NeoEpitopeFinder

| Option | Meaning |
|---|---|
| `-sample`, `-sample_id_file` | the sample, or a CSV with a SampleId column (one output file per sample) |
| `-linx_dir` | LINX's somatic output (`<sample>.linx.neoepitope.tsv`, the candidate fusions) |
| `-somatic_vcf` | PURPLE's somatic VCF (`*` is replaced by the sample id) |
| `-ref_genome`, `-ref_genome_version` | the reference FASTA (indexed) and 37 or 38 |
| `-ensembl_data_dir` | the Ensembl data cache |
| `-req_amino_acids` | amino acids either side of the novel ones (default 15) |
| `-gene_id_file` | a CSV of GeneIds to restrict the Ensembl cache to |
| `-write_trans_data` | one row per transcript (no merging of neo-epitopes with the same amino acids) |
| `-output_dir`, `-log_level`, `-log_debug`, `-threads` | as in Java (threads only split several samples) |

### NeoScorer

| Option | Meaning |
|---|---|
| `-sample`, `-sample_id_file` | the sample, or a CSV with SampleId and optional CancerType and HasRna columns (`neo_cohort.*` files with a SampleId column) |
| `-sample_data_dir` | the default of the directories below and of the output directory |
| `-neo_dir`, `-lilac_dir`, `-purple_dir` | the finder's `<sample>.neo.neo_data.tsv`, LILAC's `<sample>.lilac.tsv`, PURPLE's purity (`*` in the PURPLE directory is replaced by the sample id) |
| `-cancer_tpm_medians_file` | the cohort TPM medians per transcript ("All" and per cancer type) |
| `-cancer_type` | the sample's cancer type in that file |
| `-isofox_dir`, `-rna_somatic_vcf`, `-rna_sample` | RNA: ISOFOX's `<rna id>.isf.transcript_data.tsv` and `.isf.neoepitope.tsv`, and a somatic VCF with an RNA genotype; both directory and VCF make the sample an RNA one (TPM source SAMPLE) |
| `-cohort_trans_exp_file` | a CSV matrix of transcript TPMs by sample, in place of ISOFOX's transcript data |
| `-score_file_dir`, `-score_file_id` | the binding reference files `neo_train_<type>[_<id>].csv` |
| `-pos_weight_file`, `-rand_dist_file`, `-likelihood_rand_dist_file`, `-exp_likelihood_rand_dist_file`, `-likelihood_file`, `-flank_pos_weight_file`, `-recognition_file`, `-expression_dist_file` | the reference files one by one |
| `-check_self_recog` | an identical peptide does not count towards recognition similarity |
| `-rank_threshold` (default 0, not applied), `-sim_threshold` (default 10) | the allele-peptide rows written: likelihood rank at or below the rank threshold, or a recognition similarity at or below the similarity threshold |
| `-write_types` | `NEOEPITOPE;ALLELE_PEPTIDE` by default (`PEPTIDE` and `SAMPLE_SUMMARY` are accepted and write nothing, as in Java) |
| `-peptide_lengths` | min-max, eg `9-11` (default 8-12; Java ignores the option, NE9) |
| `-ensembl_data_dir`, `-ref_genome_version` | the transcript amino acids, to drop peptides found in an upstream transcript's protein |
| `-output_dir`, `-output_id`, `-log_level`, `-log_debug`, `-threads` | as in Java; the port uses the threads to score a sample's peptides |

Accepted and unused, as in Java: `-validation_data_file`, `-write_peptide_scores`, `-write_summary_data` (BindScorer's own
main class).

Not ported: NEO's other main classes (BindTrainer, BindScorer, the missense and utility programs under `missense/` and
`utils/`).

## Outputs

| Step | File | Content |
|---|---|---|
| finder | `<sample>.neo.neo_data.tsv` | each neo-epitope: variant, genes, transcripts, upstream, novel and downstream amino acids, NMD and coding base counts, coding bases with positions and cigars |
| scorer | `<sample>.neo.neoepitope.tsv` | each neo-epitope with its peptide count, TPM source, RNA fragments and depth, sample, expected and effective TPMs, cohort TPMs |
| scorer | `<sample>.neo.peptide_scores.tsv` | each peptide and allele: score, rank, likelihood and its rank, expression likelihood and its rank, effective TPM, recognition similarities, allele copy number and disruption |

## How it works

`finder.rs` and `epitope.rs` follow NeoSampleTask, PmNeoEpitope, SvNeoEpitope and EpitopeUtils. `score.rs` follows the
`score` package: the data loaders, NeoEpitopeData and NeoRnaData, TpmCalculator (peptides shared by the
neo-epitopes of a variant, with TPMs apportioned over shared upstream transcripts), the medians cache and the writers.
`bind.rs` holds the parts of the `bind` package that NeoScorer uses: the position weight matrices, flank weights,
random score and likelihood distributions (ranks by interpolation), binding and expression likelihoods and the BLOSUM62
recognition similarity. Each peptide's similarity to every allele's immunogenic peptides is found once and used for all
its alleles, and the peptides are scored on `-threads` threads; both leave the results as Java's.

Java's HashSet orders that reach the output (the alleles, a peptide's neo-epitope ids) are reproduced with
`hmf_common::java_hash`. The files are written when a sample is done, so a sample that fails leaves its file empty as
Java's does.

## Testing

```
CARGO_BUILD_JOBS=2 cargo test --release
```

runs NEO's Java unit tests of the ported classes (PmNeoEpitopeTest, SvNeoEpitopeTest, NeoEpitopeUtilsTest,
TpmCalcsTest, with their synthetic transcripts and reference bases) and tests of the port's own: the findings in both
modes, the binding classes, the readers and the writers' formats. The comparisons with Java are described in [REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
