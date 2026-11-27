# isofox-rs

A port of hmftools ISOFOX (Java) to Rust: RNA-seq transcript and gene expression, novel splice junctions, retained
introns and fusions from an RNA BAM, and the RNA support of NEO's fusion neo-epitopes. These are oncoanalyser's
ISOFOX_QUANTIFICATION and NEO_ANNOTATE_FUSIONS steps. The port takes ISOFOX's command line and writes the same files;
with `-java_compat` they equal Java's byte for byte. Without it, the Java faults listed in [REPORT.md](REPORT.md) are
corrected: IF1 to IF18, except IF10, IF13 and IF15, which keep Java's behaviour. IF7 (output depending on thread
timing) is corrected in both modes.

The port reproduces two versions of ISOFOX:

- **2.1** (the default): for RNA aligned by bwa-mem3, lifted over by TARS and duplicate-marked by REDUX.
- **2.0.1** (`-isofox_version 2.0.1`): for RNA aligned by STAR and duplicate-marked by GATK MarkDuplicates. Hartwig's
  README sends STAR BAMs to ISOFOX 2.0 or earlier, and upstream oncoanalyser runs 2.0.1 after STAR.

The version is chosen by, in order: the option `-isofox_version 2.0.1|2.1`; the environment variable
`ISOFOX_VERSION` (`2.0.1` or `2.1`); the BAM header (an `@PG` line of STAR, `ID:STAR` or `PN:STAR`, selects 2.0.1);
else 2.1. `isofox-rs -version` prints `Isofox version 2.1`, or `Isofox version 2.0.1` when `ISOFOX_VERSION=2.0.1`
is set. What 2.0.1 does differently is listed in [REPORT.md](REPORT.md), section 10.

## Usage

As the bioconda wrapper `isofox` takes it:

```
isofox-rs -sample <id> -functions 'TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS' -read_length 151 \
    -bam_file <rna bam> -ref_genome <fasta> -ref_genome_version 38 -ensembl_data_dir <dir> \
    -driver_gene_panel <tsv> -known_fusion_file <csv> -gene_distribution_file <csv.gz> -alt_sj_cohort_file <csv.gz> \
    -exp_counts_file <csv> -exp_gc_ratios_file <csv> -threads <n> -output_dir <dir>
isofox-rs -sample <rna id> -tumor <dna id> -functions NEO_EPITOPES -neo_dir <NEO finder dir> -bam_file <rna bam> \
    -read_length 151 -ref_genome <fasta> -ref_genome_version 38 -ensembl_data_dir <dir> -output_dir <dir>
```

JVM options (`-Xmx...`, `-Xms...`, `-XX:...`, `-D...=...`) are accepted and ignored. `isofox-rs -version` prints the
version as Java does.

Port options: `-isofox_version 2.0.1|2.1` (or `ISOFOX_VERSION`) chooses the ISOFOX version (above);
`-java_compat` (or `ISOFOX_JAVA_COMPAT=1`) reproduces every Java fault; `-java_compat_except IF5,IF12` (or
`ISOFOX_JAVA_COMPAT_EXCEPT`) all but those listed.

In 2.0.1 mode the command line is ISOFOX 2.0.1's: `-excluded_regions <file>` and `-enriched_gene_ids <ids>` are
accepted, the write types are 2.0.1's (`CHIMERIC_READ` in place of `FUSION_FRAGMENT` and `MULTI_MAP_LOCI`), and
`-single_map_qual` defaults to 255 (STAR's unique map quality). Two options of 2.1 stay accepted although 2.0.1
rejects them: `-tumor`, which oncoanalyser's NEO_ANNOTATE_FUSIONS passes (2.0.1 names NEO's file by `-sample`; the
port takes `-tumor` where given), and `-log_read_ids`. Not ported in 2.0.1 mode: the function `READ_COUNTS` and the
flags `-write_chimeric_reads` and `-write_chimeric_frags` (each ends the run with an error).

| Option | Meaning |
|---|---|
| `-sample`, `-output_dir`, `-output_id` | the sample (file prefix), the output directory, an optional extra file id |
| `-functions` | `TRANSCRIPT_COUNTS`, `ALT_SPLICE_JUNCTIONS`, `RETAINED_INTRONS`, `FUSIONS`, `STATISTICS`, `NEO_EPITOPES` separated by `;` (default the first, second and fourth) |
| `-bam_file`, `-ref_genome`, `-ref_genome_version` | the RNA BAM (indexed), the reference FASTA (indexed), 37 or 38 |
| `-ensembl_data_dir`, `-canonical_only` | the Ensembl data cache; canonical transcripts only |
| `-read_length` | the read length if known (0: taken from the reads) |
| `-exp_counts_file`, `-exp_gc_ratios_file`, `-exp_rate_frag_lengths` | expected fragment counts per gene collection and category, expected GC ratios per transcript, the fragment lengths and weights of the expected counts (default `50-1;75-1;...;550-1`) |
| `-frag_length_min_count`, `-long_frag_limit` | fragments sampled for the length distribution (default 1,000,000 with expected counts), the longest fragment (default 550) |
| `-gene_distribution_file`, `-cancer_type`, `-driver_gene_panel` | the cohort TPM medians and percentiles, the sample's cancer type in it, the driver genes (high expression reporting, novel junction genes) |
| `-alt_sj_cohort_file` | the cohort frequency of novel splice junctions |
| `-known_fusion_file`, `-fusion_cohort_file`, `-fusion_min_frags_filter`, `-fusion_skip_non_genic`, `-fusion_write_chimeric_only` | the known fusion data, a fusion cohort, the hard filter's fragments (default 2), skipping non-genic chimeric groups, writing no fusions |
| `-panel_tpm_norm_file` | panel TPM normalisation factors |
| `-specific_regions`, `-specific_chr`, `-gene_id_file` | restrict the run to regions or chromosomes; restrict NEO_EPITOPES to genes |
| `-single_map_qual`, `-gene_read_limit`, `-drop_dups` | the map quality of fragment-length and neo-epitope reads (default 60; 255 in 2.0.1 mode, where it also applies to all reads when transcript counts are not run), a per-gene-collection read limit, duplicates left out of expression (2.0.1 mode; 2.1 always leaves them out) |
| `-excluded_regions`, `-enriched_gene_ids` | 2.0.1 mode: regions whose reads are skipped; the genes counted but not otherwise processed (default RN7SL1-5 and RN7SK, which 2.0.1 also excludes with 1 kb either side) |
| `-write_type` | `FRAG_LENGTH` and `FRAG_LENGTH_BY_GENE` are written; the other write types are accepted and write nothing |
| `-neo_dir`, `-tumor` | NEO's finder output and the DNA tumor id of its file name (NEO_EPITOPES) |
| `-threads`, `-log_level`, `-log_debug`, `-log_read_ids` | as in Java; the output does not depend on the thread count |

## Outputs

`<sample>.isf.<file>`. In 2.0.1 mode the multi-mapped columns are named `LowMapQual` (gene collections) and
`LowMapQualFrags` (gene and transcript data), and the summary has `EnrichedGenePercent` before `MedianGCRatio`.

| File | Content |
|---|---|
| `gene_collection.tsv` | each gene collection (overlapping genes): fragment counts by type |
| `gene_data.tsv` | each gene: spliced and unspliced fragments, TPMs, residuals, cohort medians and percentiles, reported status |
| `transcript_data.tsv` | each transcript: fitted fragments before and after the GC adjustment, TPMs, bases covered, splice junction support |
| `summary.tsv` | QC status, total and duplicate fragments, the shares of spliced, unspliced, alternative and chimeric fragments, fragment length percentiles, median GC ratio, strandedness |
| `alt_splice_junc_unfiltered.tsv`, `alt_splice_junc.tsv` | novel splice junctions (all, and those passing: known exon deletions and duplications) |
| `retained_intron.tsv` | retained introns at exon boundaries (`retained_intron.csv`, comma-separated, in 2.0.1 mode) |
| `fusions.tsv`, `pass_fusions.tsv` | fusion candidates and those passing the filters |
| `neoepitope.tsv` | NEO_EPITOPES: fragments supporting each fusion neo-epitope and the depth at its junction ends |
| `frag_length.csv`, `frag_length_by_gene.csv` | with `-write_type FRAG_LENGTH` / `FRAG_LENGTH_BY_GENE` |

## How it works

The structure follows the Java classes (see NOTES.md for the file map). Each chromosome's gene collections are processed
on `-threads` threads: the BAM is sliced per collection (rust-htslib), reads are paired into fragments and matched against
the collection's exon regions, and fragments are counted by transcript category, novel junction, retained intron or
chimeric group. The expected counts are fitted per collection (expectation maximisation), refitted after the GC
adjustment, and TPMs worked out over all collections. Chimeric groups are completed across gene collections, then across
chromosomes, and turned into fusions.

Java writes several files from its chromosome threads and numbers fusions as threads reach them, so its output depends
on thread timing (IF7). The port buffers each chromosome's lines, runs the inter-chromosomal fusion stage in the
single-thread order, numbers fusions as one thread would and writes in that order: its output at any thread count equals
Java's at one thread. Java's hash orders that reach the output (identity-hashed exon regions, read-name and transcript-id
keyed maps) are reproduced as a JVM with `-XX:hashCode=2` iterates them (`hmf_common::java_id_hash`, `java_seq_hash`).

The 2.0.1 mode (`version.rs`) has its own fragment allocator (`allocator201.rs`: duplicates and secondary
alignments read and weighted, map-quality multi-mapping, supplementaries paired by name, excluded regions and enriched
genes), and asks `version::v201()` where the shared code differs (read classification, chimeric reads, fusion
fragments and junction types, novel junctions, the writers).

## Testing

```
CARGO_BUILD_JOBS=2 cargo test --release
```

runs ports of ISOFOX's Java unit tests (ReadUtilsTest, ReadCountsTest, TransClassificationTest, FragmentTest,
NovelJunctionsTest, GcBiasTest, MiscTest, ExpectedRatesTest's fit, FusionFragmentsTest and FusionDataTest cases) on synthetic
genes, transcripts and reads, and tests of the port's own (the findings in both modes, the hash orders, the parsers, and
each behaviour of the 2.0.1 mode in `tests201.rs`). The comparisons with Java are described in [REPORT.md](REPORT.md).
