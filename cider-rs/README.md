# cider-rs

A Rust port of hmftools CIDER 1.2 (`../cider`): the IG/TCR CDR3 sequences of a sample, assembled from the
reads at the IG/TCR loci of its BAM or CRAM, annotated with their V, D and J genes by BWA-MEM. It takes the
same options and writes the same files as the Java `CiderApplication`; with `-java_compat` every file is
byte-identical to Java's (cider.bam and the gzip files too). Faults found in the Java code are corrected
by default (findings CI1-CI11 in [REPORT.md](REPORT.md)).

## Build and run

```
cargo build --release          # target/release/cider-rs
cider-rs -sample S -bam S.bam -ref_genome_version 38 -ref_genome GRCh38.fasta -write_cider_bam \
    -threads 2 -output_dir cider/
```

As oncoanalyser calls it: `cider -Xmx… com.hartwig.hmftools.cider.CiderApplication -sample … -output_dir cider/`.
JVM options (`-Xmx…`, `-Xms…`, `-XX:…`, `-D…=…`) are dropped and the class name is accepted as the first
argument (another `com.hartwig.` class is refused, exit 2); `cider-rs -version` prints `Cider version 1.2`
and `-help` lists the options as Java does.

The BWA-MEM index image of the reference (`<ref_genome>.img`, or `-bwa_index_image`) and its `.dict` are
needed for the alignment annotation; without `-ref_genome` the alignment is skipped, as in Java. BWA is
the vendored library of hmf-common-rs (feature `bwa`, the same sources and calls as gatk-bwamem-jni); no
`-bwa_lib` is needed (the option is accepted). The IG/TCR gene table, the IMGT allele FASTA and its BWA
image, and the GRCh37 TRBJ1 patch image are compiled in (`resources/`).

## Options

| Option | Default | |
|---|---|---|
| `-sample` | | sample id, the output file prefix |
| `-bam` | required | BAM or CRAM |
| `-ref_genome` | | FASTA (CRAM decoding; with its `.dict` and `.img` the alignment annotation) |
| `-ref_genome_version` | V37 | 37 or 38 |
| `-bwa_index_image` | `<ref_genome>.img` | BWA-MEM index image |
| `-bwa_lib` | | accepted, not used |
| `-write_cider_bam` | | write `<sample>.cider.bam` with the reads used |
| `-report_match_ref_seq` | | keep sequences that match the reference |
| `-max_fragment_length` | 1000 | the margin around anchors and constant regions, the mate distance |
| `-min_base_quality` | 25 | |
| `-num_trim_bases` | 0 | bases trimmed from each read end |
| `-max_low_qual_base_fraction` | 0.1 | reads with more low quality bases are left out |
| `-max_reads_per_gene` | 100000 | downsampling limit per gene type (CI5) |
| `-primer_csv`, `-primer_mismatch_max` | , 0 | primers (tab-separated: gene, probeName, sequence, vj) to match against the sequences |
| `-threads` | 1 | worker threads (0 is taken as 1; Java fails, CI7) |
| `-output_dir` | | created if missing |
| `-log_debug`, `-log_level` | | debug logging, with `PerfStats` timings |
| `-java_compat`, `-java_compat_except CI1,...` | | Java behaviour for every finding, or all but those listed (also `CIDER_JAVA_COMPAT=1`, `CIDER_JAVA_COMPAT_EXCEPT`) |

Environment (port only): `CIDER_IMAGE_PRELOAD=0` (or 1) turns off (on) reading the BWA image through before
the alignment (on by default, see below); `CIDER_FILE_BUFFER_KB` sets the BAM read size (1024; 0 for the
htslib default).

## Outputs

In `-output_dir`:
- `<sample>.cider.vdj.tsv.gz`: the VDJ sequences with their filters, anchors, read counts, genes (V, D, J,
  supplementary), identity to IMGT (somatic hypermutation status), the layout and its support
- `<sample>.cider.locus_stats.tsv`: per locus, reads used and total, downsampled, sequences, passing sequences
- `<sample>.cider.layout.gz`: the read layouts per gene type, with each read lined up (ANSI coloured)
- `<sample>.cider.alignment_match.tsv.gz` (with alignment): the V, D, J and full alignments of each sequence
- `<sample>.cider.bam` (with `-write_cider_bam`): the matched reads, coordinate sorted as htsjdk sorts them,
  encoded as htsjdk writes them
- `<sample>.cider.primer_match.tsv` (with `-primer_csv`)

## How it works, and speed

The regions (each anchor location and constant/D region, widened by `-max_fragment_length`) are read in
Java's order; records seen before (name, pair end, reference, start) and duplicates are dropped. With 2 or
more threads one thread reads and decodes the BAM and the others screen the records (anchor location
overlap through a per-chromosome index; the BLOSUM anchor search of soft clips near constant regions and
of unmapped reads whose mate maps nearby, scored from per-codon tables). The results are taken in read
order, so every output equals that of a one-thread Java run at any thread count. Then the layouts per gene
type (in parallel), the VDJ sequences (merging per gene type pair in parallel), and the BWA-MEM alignment
on `-threads` threads.

The BWA image (5.6 GB for GRCh38) is read through once before the alignment, which then runs from the page
cache (1 s when cached, about 10 s from a USB SSD). Java maps the image the same way but does not read it:
on a Mac with little free memory its page faults made the alignment of a few hundred sequences take 80 to
180 s, in Java and in the port alike, against about 5 s after a read-through.

Whole-genome tumor / normal, 2 threads: 5.9 / 8.9 s against Java's 13.5 / 15.4 s with the image cached
(Java up to 130 s when its page faults stall), a fifth to a third of Java's CPU time, 0.6 GB peak memory
footprint against 3.6-3.9 GB (peak RSS counts the mapped image in both). BWA-MEM, the same C code as
Java's, is half of the time on the normal. See REPORT.md.

## Tests

```
cargo test --release
```

runs the ported Java unit tests (all of the CiderApplication classes: layout forest, layout builder, VDJ
builder and builder utils, annotator, read screener, BLOSUM search and similarity, read slices, alignment
blocks, somatic hypermutation) and own tests (BAM encoding, gene file, compat switches, the fast BLOSUM path
against the string path, the IMGT resources against their BWA image). The comparisons with Java are in
[REPORT.md](REPORT.md); the comparison scripts are kept outside this repository, as they name private
inputs.
