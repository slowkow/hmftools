# pave-rs

A Rust port of hmftools PAVE 2.0 (`../pave`): the transcript and protein impacts of SAGE's small
variants, with their canonical and worst impact per variant, and the annotations and filters PAVE adds
(ClinVar, the panel of normals, gnomAD frequencies, mappability, blacklists, reportability). It takes
the same options and writes the same files as the Java `PaveApplication`. With `-java_compat` the VCF
is byte-identical to Java's (hmftools master `883ab11bba`), the BGZF bytes and the tabix index too. Faults found in
the Java code are corrected by default (findings PA1-PA9 in [REPORT.md](REPORT.md); PA10 is kept as Java does).

## Build and run

```
cargo build --release          # target/release/pave-rs
pave-rs -sample TUMOR -input_vcf TUMOR.sage.somatic.vcf.gz -ref_genome ref.fasta -ref_genome_version 38 \
    -ensembl_data_dir ensembl_data -driver_gene_panel DriverGenePanel.38.tsv \
    -pon_file hmf_wgs_sage_pon_1000.38.tsv.gz -gnomad_freq_dir gnomad -clinvar_vcf clinvar.38.vcf.gz \
    -mappability_bed mappability_150.38.bed.gz -threads 8 -output_dir pave/ -output_vcf pave/TUMOR.pave.somatic.vcf.gz
```

These are the somatic options oncoanalyser passes; its germline run gives no PON or gnomAD files but
`-blacklist_bed`, `-blacklist_vcf` and `-gnomad_no_filter`. JVM options (`-Xmx...`, `-XX:...`) are
accepted and ignored, so the bioconda wrapper's command line works, and `pave-rs -version` prints
`Pave version 2.0` as oncoanalyser expects.

## Options

As Java's (`pave-rs -help` lists them):
- **Inputs:** `-sample`, `-input_vcf` (bgzipped and indexed, or plain), `-ref_genome`,
  `-ref_genome_version` (37 or 38), `-ensembl_data_dir`, `-driver_gene_panel`.
- **Annotations:** `-pon_file`, `-pon_artefact_file`, `-pon_filters` (TIER:SAMPLES:MAX_READS separated by
  ';', defaults by genome version), `-gnomad_freq_dir` or `-gnomad_freq_file`, `-gnomad_no_filter`,
  `-mappability_bed`, `-clinvar_vcf`, `-blacklist_bed`, `-blacklist_vcf`. PON and gnomAD annotations
  already in the input (its header declares PON_COUNT or GND_FREQ) are used instead of the files.
- **Sequencing type:** `-sequencing_type` ILLUMINA (default), SBX or ULTIMA: the PON filter of MSI
  INDELs found only in the technology's PON.
- **Filtering and output:** `-filter_type` PASS (default: only passing records are annotated and
  written), PANEL or ALL; `-write_pass_only`; `-set_reportable` (HOTSPOT and REPORTED flags);
  `-write_detailed` (PR_CXT); `-write_transcript_data` (`<sample>.pave.transcript.tsv`);
  `-only_canonical`.
- **Restrictions:** `-specific_chr`, `-specific_regions`.
- **Running:** `-threads`, `-log_level`, `-log_debug`, `-output_dir`, `-output_vcf`.
- **Java behaviour:** `-java_compat` gives Java's behaviour for every finding;
  `-java_compat_except PA2,...` gives it for all but those listed. The environment variables
  `PAVE_JAVA_COMPAT=1` and `PAVE_JAVA_COMPAT_EXCEPT` do the same.

## Outputs

- the annotated VCF (`-output_vcf`, else the input's name with `.pave` before `.vcf`, always bgzipped)
  with its tabix index `.tbi`;
- `pave.version` (version and build time);
- with `-write_transcript_data`, `<sample>.pave.transcript.tsv`: one line per transcript impact.

## How it works

The input VCF is read once; the records each chromosome task would process (by filter type and
specific regions) are kept, and the positions they need are collected. The resources are then read in
parallel, keeping only the entries at those positions: Java holds each chromosome's whole PON, gnomAD
and mappability data (0.8-2.7 GB of text per file), the port a few thousand entries. The chromosomes are
then processed in parallel, as Java's tasks are: genes near each variant, impacts per transcript (coding
context, protein context, HGVS strings, splice effects), the right-aligned form of INDELs with
microhomology, the phasing of variants of one local phase set, the summary impact, annotations and
filters. A variant in one of the special splice regions of the built-in resource `resources/splice/
special_splice_regions.tsv` (Java's `splice/special_splice_regions.tsv` in its jar; upstream `4adc014589`: MET exon
14, ENST00000397752, intronic bases -30 to -10 before the exon) is marked as in a splice region, as one in the
standard splice region (8 intronic and 3 exonic bases at an exon boundary) is. The records are encoded as htsjdk's VCFEncoder does and written in chromosome order to a BGZF
stream deflated by zlib at level 5 with htsjdk's block size, with the tabix index built as htsjdk
builds it, so both files equal Java's.

## Tests

`cargo test --release` runs the Java unit tests of PAVE, ported (`src/tests/`, one module per Java test
class; the tests Java ignores are ignored), and tests of the port (`src/tests/port.rs`: each corrected
finding in both modes, the resource readers, the output record). The comparisons with Java are described in
[REPORT.md](REPORT.md).

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); follows upstream master `883ab11bba` (2026-11-27) since the
review of 2026-11-28: `-java_compat` reproduces the Java built at that commit. Upstream's changes, and what was decided
about each, are in [`../UPSTREAM.md`](../UPSTREAM.md) and REPORT.md section 7.
