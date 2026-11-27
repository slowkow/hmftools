# CIDER in Rust: record of the port

A record of how `cider-rs` was checked against Java CIDER 1.2 (`hmftools/cider/target/cider-1.2-jar-with-dependencies.jar`,
the fork's HEAD), the faults found in the Java code, and timings. Options and outputs are in
[README.md](README.md).

## 1. How it was checked

Java CIDER 1.2 and the port ran on the same inputs with the same options, and every output file was
compared byte for byte (for gzip files the decompressed text too, for the BAM the
`samtools view -h` text too). Java ran with `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2` and
`-bwa_lib hmftools/libbwwwa.Darwin.dylib` (the jar holds an x86-64 library only); the port in Java mode
(`CIDER_JAVA_COMPAT=1`).

Inputs: the private whole-genome tumor/normal pair (ids TUMOR, NORMAL; 1.0 M records in the CIDER regions
of each), the public GIAB HG008 tumor's `cider.bam` from the pipeline (31,407 reads, the reads CIDER
selected), and a CRAM copy of it; GRCh38 (`GRCh38_masked_exclusions_alts_hlas.fasta` with `.dict` and the
5.6 GB `.img`), `-ref_genome_version 38`, `-write_cider_bam` as oncoanalyser runs it.

Java's identity hashes do not reach the output (its HashMaps with object keys are only looked up or
sorted), and its outputs at 1 and 2 threads were byte-identical on both whole-genome samples and on
HG008, though its read candidate list is filled in thread order (CI9). The port's outputs equal those of a
one-thread Java run at any thread count.

"identical" = every output file byte-identical: cider.bam (BGZF blocks, deflate bytes and all), the gzip
files (zlib level 6 as java.util.zip), locus_stats.

| case | options besides the inputs | result |
|---|---|---|
| wgs_tumor_t1, _t2 | whole genome, -threads 1, 2 | identical (5 files) |
| wgs_normal_t1, _t2 | whole genome, -threads 1, 2 (96 V-J joins, 89 PASS sequences) | identical (5 files) |
| hg008_t1, _t2 | HG008 cider.bam | identical (5 files) |
| hg008_noref | no -ref_genome (alignment skipped, SKIPPED_ALIGN) | identical (4 files) |
| hg008_report_ref | -report_match_ref_seq | identical |
| hg008_trim | -num_trim_bases 3 | identical |
| hg008_quals | -min_base_quality 30 -max_low_qual_base_fraction 0.2 | identical |
| hg008_frag | -max_fragment_length 500 | identical |
| hg008_primer | -primer_csv (3 synthetic primers) -primer_mismatch_max 2 | identical (with primer_match.tsv) |
| normal_down | -max_reads_per_gene 2000 (6 gene types downsampled, Kotlin Random(0)) | identical |
| normal_t0 | -threads 0 | same failure after cider.bam (identical), exit 1 (CI7) |
| normal_noref_t1 | no -ref_genome, -threads 1 | identical (3 files) |

The unit tests: every test of the Java test tree for the ported classes (LayoutForestTest,
VJReadLayoutBuilderTest, VDJSequenceBuilderTest, VdjBuilderUtilsTest, VdjAnnotatorTest,
CiderReadScreenerTest, AnchorBlosumSearcherTest, BlosumSimilarityCalcTest, ReadSliceTest, CiderUtilsTest,
SomaticHypermutationTest), and own tests; 95 in all, all passing (`cargo test --release`). Not ported:
ImgtGeneCuratorTest (the curator builds the resource files and is not part of CiderApplication) and the two
empty test classes.

### HG008: my Java run against the pipeline's output

The pipeline's HG008 CIDER output (`data/public/hg008/pipeline/cider/`, biocontainers CIDER 1.2 on the
whole HG008 BAM) and Java CIDER 1.2 of the fork run on that output's `HG008-T.cider.bam` agree: cider.bam
is byte-identical (the same 31,407 reads; CIDER selects all of them again, and htsjdk writes them the same
way), the layout file, locus_stats and alignment_match are identical, and the vdj.tsv rows are the same.
The only difference is two columns, `vImgtRefUnaligned` and `jImgtRefUnaligned`, which the fork's 1.2
writes and the pipeline's 1.2 build does not: the two "1.2" builds are not the same source. The port
writes the fork's columns. That the whole BAM and its cider.bam give the same VDJ sequences also means no
mate read outside cider.bam extended a sequence in this sample.

**In the pipeline (HG008 on Google Batch, 2026-11-27).** The port in Java mode against that pipeline
output, now on the whole HG008 BAM. It took 20 s and 0.52 GB against Java's 37 s and 6.2 GB.
- **Columns:** the vdj.tsv rows have the fork's two extra columns, as above.
- **One extra read:** cider.bam has 31,408 reads against 31,407 (flag 147 at chr14:106,356,136, MAPQ 0),
  so the IGH counts in locus_stats are one higher. Java CIDER 1.2 of the fork includes that read too, at 1
  and at 6 threads, on an IGH-locus slice of the same BAM. So the difference is the pipeline's other
  1.2 build, not the port.
- **Header:** cider.bam's `@HD VN` was the input's 1.5, while Java writes htsjdk's 1.6. htsjdk's decoded
  header does not keep the input's version. Fixed on 2026-11-27 in both modes (hmf-common-rs
  `sam_header`). On that slice the port's five outputs, cider.bam included, are then byte-identical to
  the fork's Java. Every local input said 1.6.

## 2. Findings

Each finding with a switch is corrected by default; `-java_compat` (or `-java_compat_except CIn,...`)
reproduces Java. The effect on the whole-genome pair was measured by running Java mode with one finding
corrected at a time; only CI1 and CI2 changed any output there.

**CI1: genes on ALT, patch and unplaced contigs are taken as primary-assembly locations.**
IgTcrGene.fromCommonIgTcrGene builds the gene and anchor locations with
`GenomicLocation.fromChrBaseRegionStrand(region, strand)`, which drops the assembly name of the gene file
row. Every location then counts as primary (`altAssemblyName == null`), so CiderGeneDataLoader's filter of
non-primary locations does nothing. 64 rows of igtcr_gene.38.tsv have `primaryAssembly` false; their
coordinates are positions within an ALT, patch or unplaced contig, written against the primary
chromosome's name (for example TRBV genes of the HSCHR7_2_CTG6 ALT at chr7:0-1 Mb). Java adds 40 anchor
locations at those positions (466 instead of 426; 526 query regions instead of 485) and takes reads there,
in unrelated genes on chr7, chr2 and chr16 below 3 Mb, as V or J anchor reads by alignment alone. Effect on
the whole-genome pair: 2,076 (normal) and 1,437 (tumor) such reads in cider.bam, the locus_stats read
totals inflated (normal TRB 8,292 against 6,536; IGH and IGK by about 1.5%), 90 and 78 extra layouts, and
5 and 4 spurious TRB sequences, all filtered NO_VDJ_ALIGNMENT; the PASS sequences are unchanged. For
GRCh37 the patch genes TRBJ1-1 to 1-6 (on 7_gl582971_fix) land at chr7:~990 kb, a gene-dense region.
Corrected: the assembly name is kept, as `toCommonIgTcrGene` and the curator write it.

**CI2: the adapter trim takes the read's 3' end as its 5' end.** CiderUtils.calculateAdapterSequenceTrim
calls `getFivePrimeUnclippedPosition(alignmentStart, cigar, read.readNegativeStrandFlag)` whose third
argument is `forwardStrand`, so for either strand it gets the read's 3' end, and the insert length tested
against the read length is that of the mate's 5' end to the read's 3' end. Corrected: the read's 5' end.
Effect: one read (of a 71-base fragment, in a CI1 region) in the normal's cider.bam; none in the tumor.

**CI3: a V anchor ending at the last base of a sequence is rejected.** AnchorBlosumSearcher.tryMatchWithBlosum
tests `anchorEnd >= dnaSeq.length` to reject a V anchor that would run past the sequence (the comment: "we
don't want to allow removing the conserved C"); an anchor that ends exactly at the end fits. Corrected: `>`.
Only unmapped-read searches reach the end of a sequence; no effect on these samples.

**CI4: the closest alignment block is measured from the block's start.** In
CiderReadScreener.extrapolateReadOffsetAtRefPosition the distance to a block that ends before the position
is `referencePos - blockReferenceStart`, not `- blockReferenceEnd`, so a position in a deletion (or intron)
is extrapolated from the following block though it lies next to the one before (60M10D90M: 3 bases into the
deletion gives read offset 52, not 62). Corrected: from the block's end. No effect on these samples.

**CI5: downsampling at exactly the limit.** VJReadLayoutBuilder.buildLayouts downsamples when the read count
is `>= maxReadCountPerGene`; at the limit the fraction kept is 1, so no read is dropped but locus_stats
reports `downSampled` true; the option says reads are downsampled when they exceed the limit. Corrected: `>`.

**CI6: mergeAnchors keeps the first anchor.** VdjBuilderUtils.mergeAnchors means to keep the read-matched
anchor with more reads but calls `maxOf(anchor1, anchor1, comparingInt(numReads))`. Corrected:
`maxOf(anchor1, anchor2, ...)` (the first on a tie). No effect on these samples: the merge lists are sorted
by read count, so the first anchor nearly always has the most reads.

**CI7: `-threads 0` fails.** buildLayouts creates `Executors.newFixedThreadPool(threadCount)`, which throws for
0 or less, after cider.bam is written; the option says "0 or 1 not multi-threaded". Corrected: taken as 1.
Java mode writes cider.bam and exits 1 as Java does.

**CI8: query regions with the same start are merged.** The regions go into a TreeSet whose
GenomeRegion.compareTo compares the chromosome (by rank) and the start only, so of two regions with the
same start the one added later is dropped even when it reaches further. Corrected: such a region is kept
(only identical regions are merged). No region was dropped for GRCh38.

**CI9 (no switch): the read candidates are collected in thread order.** With `-threads` above 1 the
records are screened by consumer threads appending to synchronized lists, so the order of the read
candidates, of the unmatched reads and so of the mate candidates depends on timing. It reaches the output
only through ties (two records of the same read and pair end with the same aligned position, quality and
end; mates equally close to two anchors), and Java's 1- and 2-thread outputs were identical here. The port
always takes the one-thread order (cider.bam is coordinate-sorted by htsjdk's comparator, which does not
depend on it).

**CI10 (no switch, Java kept): a sort whose result is dropped.** LayoutForest.checkAndSealNodes calls
`readsToReassign.sortedBy { alignedPosition }` and discards the result, so the reads go back in the order
collected. Sorting them ascending would add reads against the order the forest expects (largest aligned
position first), so the intent is unclear; Java's behaviour is kept.

**CI11: a stale consensus.** ReadLayout clears its consensus cache only when a likely base changed, so a
read that lengthens or shifts the layout with N bases only leaves the previous consensus (shorter, or
shifted) in place. Only the mate extension adds reads after the consensus was taken. Corrected: cleared
after every read added. No effect on these samples.

Noted, not changed:
- A record that makes a consumer thread throw (a read with a base other than ACGTN in a BLOSUM search
  throws IllegalArgumentException; an unmapped read with a mapped position but an unmapped mate fails a
  `require`) kills that thread; with one thread the reader then waits for ever on its full queue (no
  timeout). The port stops with the error.
- The mate's end is inferred as its start plus this read's length (isUnamppedReadRelevantToAnchorLoc,
  ...ConstantRegion), though the MC tag is at hand; a documented approximation.
- findSequenceOverlap needs more than `minOverlappedBases` (21 for 20) matching bases; the intent is unclear.
- `-primer_csv` is read tab-separated (PrimerTsvFile).
- The GATK BwaMemAligner setters of the match and mismatch scores leave the scoring matrix as
  mem_opt_init filled it; CIDER's 1/-4 are the defaults, so the matrix is consistent.

## 3. Timings and memory

`/usr/bin/time -l` on the same Mac (M-series, 20 cores, 128 GB, heavily shared: about 14 GB of swap in
use), inputs and the BWA image on a USB SSD. Each case ran the port (Java mode) and then Java at once
so that Java sees the page cache the port leaves; all outputs were identical.
CPU = user + sys; footprint = macOS peak memory footprint (peak RSS also counts the pages of the mapped BWA
image: Java 6.7-9.5 GB, the port 6.1-6.4 GB on the whole genome).

| case | threads | Java real / CPU / footprint | port real / CPU / footprint | speed-up |
|---|---|---|---|---|
| HG008 cider.bam | 1 | 3.9 s / 10.9 s / 0.94 GB | 1.9 s / 1.9 s / 0.13 GB | 2.1x |
| HG008 cider.bam | 2 | 3.4 s / 10.8 s / 0.76 GB | 1.6 s / 1.9 s / 0.14 GB | 2.2x |
| WGS tumor | 1 | 17.9 s / 38.4 s / 3.92 GB | 9.0 s / 8.9 s / 0.61 GB | 2.0x |
| WGS tumor | 2 | 13.5 s / 38.4 s / 3.69 GB | 5.9 s / 8.9 s / 0.60 GB | 2.3x |
| WGS normal | 1 | 22.9 s / 42.6 s / 3.96 GB | 14.9 s / 14.9 s / 0.62 GB | 1.5x |
| WGS normal | 2 | 15.4 s / 42.0 s / 3.64 GB | 8.9 s / 14.8 s / 0.62 GB | 1.7x |

Where the time goes (2 threads, tumor / normal; the port's `-log_debug` PerfStats, Java's log times):

| step | Java | port |
|---|---|---|
| reading the BAM (1.0 M records in 526 regions, 752 k screened) | 6.9 / 6.0 s | 1.8 / 1.9 s |
| cider.bam | (in the above) | 0.13 / 0.14 s |
| layouts | 1.7 / 1.8 s | 0.53 / 0.55 s |
| VDJ sequences | 1.4 / 1.3 s | 0.14 / 0.13 s |
| alignment (BWA-MEM, genome and IMGT) | 2.9 / 5.7 s | 3.2 / 6.1 s (0.7 s of it reading the image) |
| JVM start, gene data | about 0.5 s | 0.01 s |

BWA-MEM is the same C code in both, so the alignment cannot get faster than Java's; it is about half of the
port's time on the normal (more sequences to align), and scales with threads (5.7 / 11.7 s at 1 thread).
The rest is 3.5-4 times faster than Java at the same thread count. The BAM pass: at 1 thread 2.2 s (a
`sample` profile: about a quarter disk reads, 15% inflate, the rest record decoding, conversion and
screening, after the anchor locations were indexed by chromosome and position: the per-read scan of all 466
anchors had been 40% of the pass); at 2 threads one thread reads and decodes while the other screens.

**The BWA image and the page cache.** Java's library and the port map the 5.6 GB image the same way
(`mmap`, MAP_SHARED). On this machine page faults into the mapping made the alignment of 271 sequences take
80 to 180 s, in the port as in Java (Java took 130 s for the tumor right after a port run; 14 s for the
same input minutes earlier), even when a sequential read of the file a minute before had found it cached.
The port therefore reads the image through before aligning (on by default; `CIDER_IMAGE_PRELOAD=0` turns
it off): 0.6-1 s when cached, about 10 s when not, and then the alignment takes 3-6 s. The shared module
does this on FUSE and NFS file systems for every user; on Google Batch through a bucket mount, where the
pipeline's Java CIDER took 55 min at 3% CPU, each page of the image is fetched by its own read, which this
avoids. The comparisons above were made with the image cached and Java's pages resident (Java's best case).

Memory: the port holds the screened reads (about 0.7 M unmatched mapped reads are kept for the mate
search, as Java does) in 0.6 GB; Java's footprint is 3.6-4.1 GB (its heap; `-Xmx8g`).

## 4. Left out

- The BLAST code (AlignmentUtil.runBlastn, blastnMatchtoGenomicLocation): CIDER 1.2 aligns with BWA-MEM
  only and has no BLAST option.
- The IMGT gene curator (`curator/`), a separate tool that builds the resource files the port compiles in.
- `-bwa_lib` is accepted and not used (BWA is linked in).

## 5. Open points

- GRCh37: the 37 gene file, IMGT image and the TRBJ1 patch alignment are ported and their resources load
  (unit test), but there was no GRCh37 input to compare with Java.
- CRAM input: all outputs but cider.bam identical on a CRAM 3.0 copy of the HG008 cider.bam; cider.bam has
  the same records in the same order with the tags in another order (htsjdk's CRAM decoder puts RG after the
  first tag, htslib last; htslib is told not to add MD/NM, which htsjdk does not make up, and the tags are
  re-encoded and the bin computed as htsjdk does for records not read from a BAM). Java 1.2 (htsjdk 2.23)
  cannot read CRAM 3.1. On that sparse CRAM, htslib's reference checksum of each slice's span on every
  region query takes 46 s (Java 60 s); a whole-genome CRAM has dense slices.
- The tumor and HG008 join no V and J layouts; the normal (96 joins) and the unit tests cover the joining.
