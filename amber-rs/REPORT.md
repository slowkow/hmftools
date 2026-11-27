# AMBER in Rust: record of the port

A record of how `amber-rs` was checked against Java AMBER 4.4 (this repository, `mvn -pl amber -am
package`), the faults found in the Java code, and timings. Options and outputs are in
[README.md](README.md). The segmentation of the BAFs is the shared crate's (`hmf-common-rs`, its own
REPORT.md, findings SG1-SG3).

## 1. How it was checked

Java AMBER and the port ran on the same inputs with the same options, both at 2 threads unless stated.
The port ran in Java mode (`-java_compat`). Every output file was compared after decompression, and in
every case the compressed gzip and BGZF bytes and the tabix index were identical too. `amber.version`
was compared by its version line, as it holds the build time.

The core's outputs do not depend on hash order or thread count (section 3), so the Java references used
the default JVM. The tumor-only analysis' Java references used `-XX:+UnlockExperimentalVMOptions
-XX:hashCode=2` (finding AM54).

### Inputs

- A private whole-genome tumor/normal pair (about 88x and 115x; BAMs of 96 and 120 GB).
- A small private tumor/normal pair: chr21, chr22 and parts of chr1-10.
- The public GIAB HG008 SBX and Ultima chr21 slices, as BAM and as CRAM.
- Synthetic inputs: target BEDs, a VCF loci file, a contamination mix, and tumor-only data sets of up to
  6.2 million sites.

### Whole genome

| Case | Java time, peak memory | Rust time, peak memory | Result |
|---|---|---|---|
| Tumor/normal, `-write_tumor_data` | 51.6 min, 9.42 GB | 12.6 min, 1.62 GB | identical, every file |
| Germline-only | 33.9 min, 8.84 GB | 7.3 min, 1.28 GB | identical |
| Tumor-only, `-write_tumor_data` | 29.2 min, 8.47 GB | 5.5 min, 3.16 GB | identical |
| Tumor-only, `-tumor_only_min_support 0` | 45.0 min, 8.94 GB | 6.0 min, 3.28 GB | identical |

The tumor/normal case had:
- 1.67 million heterozygous and 3.41 million homozygous sites;
- 113 regions of homozygosity and 31,542 contamination sites;
- the .pcf, and the VCF with its index.

The times are those of the merged crate, at 2 threads, with the BAMs on a USB disk. Earlier runs of
the core alone, made while a Java run read the same disk, took up to 23 minutes for the tumor/normal
case, so disk contention can double the time.

### Small inputs and options

All of these cases are identical to Java. "TN" is the small pair in tumor/normal mode.

| Case | Java | Rust |
|---|---|---|
| TN, `-write_tumor_data` | 93.7 s, 5.64 GB | 18.8 s, 0.60 GB |
| TN, 1 thread (Java at 1 and 2 threads identical, Rust too) | 114.7 s, 5.70 GB | 36.5 s, 0.60 GB |
| Germline-only | 36.3 s, 5.03 GB | 10.9 s, 0.60 GB |
| Tumor-only, and at min support 0 | 45.2 s, 5.40 GB; 77.6 s, 5.99 GB | 12.9 s and 16.3 s, 3.13 GB |
| `-write_unfiltered_germline`, germline-only and TN, one chromosome each | 20.6 s; 35.4 s | 6.0 s; 9.5 s |
| `-specific_regions` (three regions) | 19.9 s | 6.4 s |
| `-target_regions_bed`, synthetic, one chr21 line out of order (AM5) | 4.9 s | 1.9 s |
| Changed thresholds: base and map quality, depth and het AF limits, tumor min depth, `-position_gap 500` | 125.1 s | 38.7 s |
| Two reference BAMs | 78.9 s | 23.8 s |
| Loci as a VCF (FILTERed records, SNPCHECK flags) | 54.3 s | 15.7 s |
| HG008 SBX, BAM and CRAM, `-sequencing_type SBX` | 35.1 s; 75.8 s | 16.8 s; 22.6 s |
| HG008 Ultima, BAM and CRAM, `-sequencing_type ULTIMA` (UQ filter at 14,844 of 17,298 tumor sites) | 98.5 s; 222.3 s | 23.6 s; 64.0 s |
| Contamination: the small tumor with 5% of the HG008 tumor's chr21 reads (contamination 0.055; Java identical over two runs) | 43.6 s | 17.7 s |
| `-threads 0` (AM7): the same files written, the same failure | 58.5 s | 16.8 s |

The default mode gives the same output as Java on every real input: no correction changes the small
pair, the SBX slice or the whole-genome normal. Only synthetic cases change (AM5, and the tumor-only
cases below).

### Tumor-only analysis

The purity analysis and LOH identification were also checked step by step. A Java harness, compiled
against the AMBER jar, runs LohIdentifier and the steps of TumorOnlyPurityAnalysis on a
`-write_tumor_data` file and writes every intermediate value:
- every site's LOH state;
- every search level's counts and score;
- every peak's bands, gnomAD means, AUCs and class;
- the BAF selection.

A Rust driver writes the same.

| Input | Result |
|---|---|
| 9 synthetic cases at min support 2 and 0: LOH from a chromosome start, in the middle, to the end, across a centromere gap and through HLA; contamination 3% and 10%; purity 0.15; low copy-number peaks; 48 sites; 1500x; sites inside centromere gaps, on X and Y, two alt alleles at a position | 54 of 54 files identical to Java |
| The small tumor at min support 0 and 2 (179,254 and 49,271 points) | identical; the harness also reproduces Java AMBER's own files |
| Whole-genome-sized synthetic data: 6.2 million sites, 10% contamination | identical, including gnomAD means summed over 885,000 points in Java's HashSet order |

Java's harness took 1,542 s and 6.6 GB on the whole-genome-sized input; the Rust driver took 8.2 s. The
purity analysis itself takes 1.3-1.6 s for 3-6 million points.

### Library checks

Each check uses a small Java program compiled against the jar.
- **commons-math3 PoissonDistribution:** 756,800 probabilities and cumulative probabilities, the ones
  AMBER's contamination model and homozygosity test use, are bit-identical.
- **TumorContaminationModel:** every level of 3,000 synthetic allele-count maps is bit-identical.
- **BinomialDistribution.cumulativeProbability:** all 5,111,983 values of the peak search's grid are
  bit-identical.
  - The match needs FastMath's own `exp` and `log`, ported operation by operation with the library's
    tables into `hmf_common::java_fastmath`.
  - With correctly rounded `exp` and `log` instead, 89 binomial values and 343 Poisson values differ in
    the last bits.
- **HashSet iteration order** of PositionEvidence keys, including tree bins and resizes: identical on 6
  key sets of 20 to 300,000 keys.
- **DecimalFormat("0.0000"):** identical to Java on 4,006,000 BAF-like values.

### Unit tests

`cargo test --release`: 115 pass, and 5 are ignored because they need files written by Java.

Ported Java tests:
- the core: PositionEvidenceTest (CIGAR walk, fragment sync), TumorEvidenceTest,
  BaseDepthIntersectFilterTest, CosanguinityAnalyserTest, RegionOfHomozygosityFinderTest,
  TumorContaminationModelTest;
- all of the purity package's tests;
- the end-to-end scenarios of RunAmberTest: twoChromosomes, newSegmenter, targetedMode and
  handleEmptyData. They run on the Java test's synthetic BAMs.

blacklistedPointsInTargetedMode cannot be ported, as AMBER 4.4 has no such option. The port's own tests
cover each correction, the LOH segments (at either chromosome end, across or partly over a centromere
gap) and the output files.

## 2. Findings

Each is corrected by default and reproduced with `-java_compat`.

### The core (AM1-AM7)

- **AM1. In genome 38 the homozygosity exclusion BED never applies.**
  - *Java:* RegionOfHomozygosityFinder compares the chromosome name without its prefix ("1") with the
    BED's names ("chr1"). So none of the BED regions is excluded: chr1 and chr9 heterochromatin, the
    chr16 gap, and the short arms of 13, 14, 15, 21 and 22. Only the centromere regions apply.
  - *Corrected:* names are compared without the prefix.
  - *Effect:* 694 sites of the 38 loci file lie in those regions outside the centromeres. The
    whole-genome normal has no homozygous region there.
- **AM2. The CIGAR walk counts hard clips as read bases and mishandles a leading insertion.**
  - *Java:* ReadCigarState adds the length of every element that does not consume the reference to the
    read index, hard clips and padding included. It also starts the walk at element 0, even when that
    element is an insertion.
  - *Effect of the fault:* a primary read with a leading hard clip gives the base H positions further
    along. A leading insertion gives the inserted base at the alignment start, and a deletion after it
    is missed.
  - *Corrected:* elements without read bases do not move the index, and the walk starts at the first
    element that consumes the reference.
  - *Effect on the data:* 1 read of 54.7 million in the small tumor, none in the SBX or Ultima slices.
    It matters for aligners that hard-clip primary alignments.
- **AM3. A read base outside ACGTN ends the run.**
  - *Java:* PositionEvidence.equalsRef and equalsAlt use AmberBase.valueOf, which throws for IUPAC codes
    and '='.
  - *Corrected:* such a base supports neither allele.
  - *Effect:* the SBX slice has 6 reads with R or M bases, but none at a site with enough base quality.
    On whole-genome SBX data Java would fail once such a read meets a site.
- **AM4. A homozygous region starting at a chromosome's first site is never reported.**
  - *Java:* createRegionIfPassFilter rejects a start index of 0 (`<= 0`), but -1 is the "no region"
    marker.
  - *Corrected:* `< 0`.
- **AM5. A target BED whose lines for a chromosome are not together keeps only the last group.**
  - *Java:* BedFileReader.loadBedFileChrMap starts a new list at each change of chromosome and puts it
    into the map over the earlier one.
  - *Corrected:* the groups are appended.
  - *Effect:* in a synthetic BED with one chr21 line after chr22, Java keeps 35 chr21 BAF sites;
    corrected, 1,131.
- **AM6. An invalid configuration exits with status 0.**
  - *Java:* main ignores run()'s result, so for example ids and BAMs of different counts log an error,
    write no output, and exit 0.
  - *Corrected:* exit 1.
- **AM7. `-threads 0` fails after the BAF file.**
  - *Java:* the option's description says 0 means single-threaded, but `newFixedThreadPool(0)` throws.
    The run ends without the .pcf and the homozygous-region file.
  - *Corrected:* 0 is taken as 1. COBALT has the same fault (CO5).

### Tumor-only purity analysis and LOH (AM51-AM55)

- **AM51. The LOH moving window is off by one at the chromosome start.**
  - *Java:* LohIdentifier sets the rates of the window centre only once 101 sites have been added and one
    dropped. The first full window is never used, so site 50 never gets rates. And since a segment is
    widened by 50 sites, the first valid site of a chromosome can never be recovered. The last window is
    used, so the two ends of a chromosome are treated differently.
  - *Corrected:* the centre is set from the first full window.
  - *Effect:* for an LOH segment starting within the first 51 sites of a chromosome, the every-5th
    recovered sites shift by one (the same density, different sites). On synthetic data with LOH from
    two chromosome starts, about 1,180 BAF lines change.
- **AM52. A centromere gap partly inside an LOH run is subtracted whole.**
  - *Java:* the full length of every overlapping gap is subtracted from the run's length.
  - *Corrected:* only the overlap is subtracted.
  - *Effect:* none with the standard loci file, which has no site inside a centromere gap.
- **AM53. The AUC overflows an int.**
  - *Java:* CategoryEvidenceIntegral multiplies two int counts, which wraps once their product passes
    2^31. That is the size of whole-genome categories, for example 50,000 peak points times 43,000
    baseline points.
  - *Corrected:* the product is taken in double.
  - *Effect:* on whole-genome-sized synthetic data with 10% contamination, Java's mutation AUC comes out
    0.823 instead of 0.999, just under the 0.825 bound, and Java reports no contamination. Corrected,
    the contamination is 0.101.
  - This is the most consequential finding: a contaminated tumor-only sample can be reported as clean.
- **AM54. The peak's gnomAD means depend on identity hashes.**
  - *Java:* the frequencies are summed in the iteration order of a HashSet of PositionEvidence, whose
    hash includes enums hashed by identity. So the last bits of the means can differ from one JVM run to
    the next. The means are compared with the baseline ±0.2, so the outcome changes only at that
    boundary.
  - *Java mode:* reproduces the order under `-XX:hashCode=2`.
  - *Corrected:* sums in genome order.
- **AM55. gnomAD frequencies are looked up by position only.**
  - *Java:* at a position with two sites (two alt alleles), the binary search returns either.
  - *Corrected:* the site with the evidence's ref and alt.
  - *Effect:* none with the standard loci file, which has one site per position.

### Kept as Java

These are not faults, or their intent is unclear.
- The homozygous-region file and the QC's UniparentalDisomy write chromosomes without the prefix in 38.
- `-write_unfiltered_germline` in tumor/normal mode makes every site a heterozygous site for the tumor
  analysis.
- An N (reference skip) in the CIGAR moves the walk like a D, but is not a deletion. The base before the
  skip is counted at sites inside it; this applies to RNA only.
- Contamination at a homozygous site needs a normal ref support above 7. The constant is named as a
  minimum normal depth.
- A TSV loci file without a Frequency column gives frequency 0, but a VCF gives 0.5. Only the tumor-only
  analysis reads it.
- The `-write_tumor_data` file gets the eight LOH columns only when its first row has LOH data. So in
  tumor-only mode it can have rows of 11 and of 19 columns.
- The SNV-type classifier takes the alt as the minor allele up to a VAF of 0.65, not 0.5.
- Log messages: the gnomAD check's trace log swaps "high" and "low", and one message reports the
  captured count where it means the points of sufficient depth.

## 3. Hash order and thread count

No output of the core depends on hash iteration order or thread timing. Java at 1 and 2 threads, and
Java twice on the contamination case, give identical files.
- Every output list is sorted by position before it is written.
- BaseDepthIntersectFilter's HashSet keeps Java's order (its hash has no identity hash).
- The contamination model's sum follows Java's bucket order. Only keys sharing a bucket can come out in
  another order, and the level was bit-identical in every check.
- A region task counts only its own sites, so counts do not depend on the thread count.

The tumor-only analysis depends on identity hashes in one place, AM54.

## 4. Timings and memory

At the same 2 threads, Rust is 3-5 times faster than Java on the small cases. On the whole genome it is
4.1-7.5 times faster:
- tumor/normal: 12.6 against 51.6 min;
- germline-only: 7.3 against 33.9 min;
- tumor-only: 5.5 against 29.2 min, and at min support 0, 6.0 against 45.0 min.

Rust uses 2-39 times less memory: 0.6 GB on the small pair and 1.3-1.6 GB on the whole genome, against
4-9 GB for Java. The tumor-only memory, 3-4 GB, is mostly the 6.3 million loci sites held with a string
per field.

The whole-genome runs were bound by the USB disk holding the BAMs: twice as many profile samples in
`read` as in decompression. Each BAM reader asks htslib for 1 MB blocks; `AMBER_FILE_BLOCK_KB` changes
it, and 0 gives htslib's default.

A small `-position_gap` is slow in Java (chr21 at 500: 125 s). Each region task is an index query, so
close tasks decompress the same reads again. The default gaps of 4,000 (BAM) and 10,000 (CRAM) give long
tasks.

### Reading through a bucket mount (2026-11-27)

On Google Batch, with the BAMs read through gcsfuse, the first version took 10.8 min on HG008 at 13% CPU
use. COBALT and BamTools, which read sequentially, took 3-4 min. Java AMBER took 12.2 min.
- **Cause:** Java reads each region task with its own index query: 22,829 tasks on the germline sites,
  covering 2,672 Mb. With 16 threads taking tasks in turn, each reader jumps through the file, and the
  bucket mount serves each query as a new ranged request; gcsfuse also shortens the requests of a reader
  that seeks.
- **Change:** a queue item is now a run of consecutive tasks of one chromosome, up to 10 Mb. A BAM's item
  is read with one index query over its span, and each read is given to the tasks whose regions it
  overlaps. The test is htslib's own for a query `[start - 1, end)`: read start < end, and `bam_endpos` >
  start - 1.
- **Why the counts cannot change:** a query returns the BAM's reads passing that test, in file order. One
  query over a run of tasks therefore returns a superset of each task's reads, in the same order, and the
  dispatch picks out exactly the reads a task's own query would return.
- **Gaps and CRAM:** a gap above 256 kb between tasks starts a new query. Gaps up to that size are about
  160 Mb of the genome (5%), read and discarded. A CRAM keeps one query per task, because htslib's CRAM
  iterator filters with CRAM's own alignment end, but one reader still takes the consecutive tasks.
- **Settings:** `AMBER_STREAM_CHUNK_MB` and `AMBER_STREAM_GAP_KB` set the two sizes, and
  `AMBER_STREAM_GAP_KB=0` gives one query per task.
- **Checks:** all 17 small cases are byte-identical to Java with the change, including the CRAM cases.
  So are three more runs with 1 Mb items and with 5 kb and 20 kb gap splits.
- **Speed:** on local files the two modes take the same time (small cases back to back: 24.0/23.7 s
  against 25.9/24.2 s, and 28.1/25.3 s against 27.5/26.3 s). With `-position_gap 500` the change is
  faster, 10.7 s against 39 s, as close regions are no longer decompressed again.

## 5. Left out, and open points

- `-bam_validation` is accepted and not applied, as htslib has no record validation like htsjdk's.
- A CRAM needs `-ref_genome`. Without it the port stops with an error, as Java fails too, rather than
  let htslib download the reference from the EBI.
- Java's exact log text (timestamps, `Double.toString`, stack traces) is not reproduced; the messages
  and their order are.
- FastMath.log1p is not ported. The correctly rounded log1p is used where commons-math3 takes it (logGamma
  below 8), and every value checked is bit-identical anyway.
- Tumor-only memory could drop several-fold by interning chromosome names and storing bases as enums in
  the loci sites.
