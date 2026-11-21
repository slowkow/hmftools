# REDUX in Rust: record of the port

A running record, kept for a later announcement and for reports to the Java maintainers. Status and
measurements are in [PLAN.md](PLAN.md).

## 1. How it was checked

Java REDUX 2.1 (this repository, `mvn -pl redux -am package`) and the port were run on the same input
with the same options. The BAMs were compared record by record with `redux-rs compare_bams`; the other
files with `cmp`.

Inputs:
- **Illumina, private:** a slice of a whole-genome tumor BAM turned back into FASTQ and realigned with
  `bwa mem -Y`, as the pipeline does, so that the input is pre-REDUX (25.8 million reads). Also split
  into two "lane" BAMs.
- **UMIs, simulated:** the same reads with UMIs added to the read names by `simulate_umis.py`: families
  with one- and two-base errors, split families, duplex halves swapped between strands, poly-G tails,
  variable lengths. No real UMI data was available.
- **SBX and Ultima, public:** chr21 slices of GIAB HG008 tumor and normal, pre-REDUX.

| Input | Run | Records | Result |
|---|---|---|---|
| Illumina | 1 thread | 27,684,550 | BAM and all four side files identical |
| Illumina | two lane BAMs | 27,684,550 | identical |
| Illumina | 8 threads | 27,684,550 | identical apart from the order of 1,004 unplaced reads (section 2) |
| Illumina, chr21 + chr22 + unplaced reads realigned | 16 threads | 114,593,080 | identical apart from the order of the unplaced reads; 40 million reads unmapped in the first pass |
| Illumina whole genome, vendor (DRAGEN) alignments, 86 GB | 1 thread | 1,780,126,609 | BAM and all four side files identical |
| The same | 16 threads | 1,780,126,609 | identical up to record 1,706,028,352; then about 500 duplicate flags differ on one very deep decoy contig, where Java depends on thread timing (section 2), and with them `duplicate_freq.tsv` |
| The same, on Linux (x86-64) | 16 threads | 1,780,126,609 | the same result as on macOS: identical up to the same record, then duplicate flags on the deep contig; BQR, jitter and microsatellite files identical. Rust 26 min and 1.6 GB, Java 50 min and 35.7 GB, local SSD |
| Two smaller vendor inputs on Linux | 4 and 16 threads | 114 M in the larger | BAMs identical apart from the order of unplaced reads; Rust 2.5 and 3.2 times faster |
| Illumina, 6.6 M reads | 14 option runs: no consensus, `-drop_duplicates` (with and without consensus), `-skip_duplicate_marking`, `-unmap_mt`, `-unmap_alt_decoy`, `-std_chr_only`, `-skip_unmapping`, `-skip_fully_unmapped`, `-specific_regions`, `-specific_chr`, `-read_length`, `-output_id` | 1.6 to 7.0 M | all identical |
| Illumina + UMIs | `-umi_enabled`; duplex with `+` and `_` delimiters; `-umi_variable_length`; `-umi_defined_ids`; `-umi_base_diff_stats` | 27.2 M | BAM, side files and UMI statistics files identical; 586,632 dual-strand groups in the duplex run |
| SBX tumor | 1 and 8 threads | 23,819,453 | identical |
| SBX normal | 1 and 8 threads | 11,013,307 | identical |
| SBX tumor | no consensus; `-sbx_max_dup_dist 0`; a region | 5.0 to 23.5 M | identical |
| Ultima tumor | no consensus, 1 and 8 threads | 21,226,640 | identical, with `ultima_stats.tsv` |
| Ultima normal | no consensus | 39,774,627 | identical |
| `-bqr_jitter_msi_only` | Illumina, SBX, Ultima; 1 to 64 threads; with MSI model files; `-write_msi_site_file` | | `bqr.tsv`, `jitter_params.tsv`, `ms_table.tsv.gz`, MSI prediction and site file identical |

Unit tests: the Java tests of the ported classes are ported (184 Rust tests in all).

## 2. Where Java's output depends on the thread count

None of these is a difference between the port and Java at one thread.

- **Partitions.** The genome is split into `2 x threads` partitions. Duplicates are not matched across a
  boundary, and the BQR block grid starts at each partition's first read, so `bqr.tsv` and a few
  duplicate groups change with the thread count. The port splits the same way and matches Java at the
  same thread count.
- **Unplaced reads.** With several threads all workers append to one synchronised writer, so the order
  of the fully unmapped reads at the end of the BAM depends on timing. The port writes Java's
  one-thread order.
- **`duplicate_freq.tsv` labels.** For group sizes above 10 the size written is the first one seen for
  its rounded bucket. Statistics are kept per thread and merged, so the label depends on which thread
  took which partition. The port gives the one-thread label.
- **`umi_coord_freq.tsv`.** `UmiStatistics.merge` keeps the first thread's maximum for a row that exists
  in both, instead of the larger, so 14 of 2,038 rows differed between Java at 8 threads and at 1.
- **`ultima_stats.tsv` row order** is HashMap order and, with several threads, depends on the same.
- **Duplicate groups in very deep regions.** Each thread keeps one `ReadCache` for all the partitions it
  takes, and the cache's size-check counter (`mCheckSizeReadCount`, one check every 10,000 reads) is not
  reset between partitions. When the cache holds 100,000 reads or more, the check decides when half of
  it is purged early, and reads that arrive later for a purged position form a new group. So the groups
  in a very deep region depend on how many reads the same thread cached before it: with one thread, on
  every earlier partition; with several threads, on which partitions that thread happened to take. Seen
  on a whole genome, on a 33 kb decoy contig holding 7.5 million reads: adding 31,000 unrelated reads on
  chr1 to the input changed the duplicate groups on the decoy at one thread. The port gives Java's
  one-thread result at any thread count: a partition whose cache nears the threshold waits for the read
  counts of the earlier partitions.
- **Interim BAM of unmapped reads.** Written by several threads in timing order, then sorted by
  position; reads at the same position and strand keep the timing order.

## 3. Java behaviour worth reporting

### Failures
- **Ultima with `-form_consensus`:** `BaseBuilder.determineBaseAndQual` returns an invalid base for
  Ultima ("TODO: decide for Ultima") wherever the reads of a duplicate group disagree, and the BAM
  encoder then throws "Bad base passed to charToCompressedBaseLow". Two reads at one position that
  differ in one base reproduce it. The port stops at the same read with a message.
- **`-umi_base_diff_stats`:** `UmiGroupBuilder.captureUmiGroupStats` casts a `ReadInfo` to `SAMRecord`,
  so the first fragment without a duplicate ends the run with a ClassCastException. The option cannot
  complete on ordinary data.
- **Duplex UMI without the delimiter:** `splitUmi` throws. A delimiter longer than one character is
  mishandled (the second half starts one character after the delimiter's start).
- **Ultima with `-ultima_drop_qual_tags` and BQR on:** null pointer in
  `UltimaBamUtils.extractT0Values`, since BQR needs the tag that was dropped.
- **Ultima read without `t0`, or with a short one:** null pointer or index error in `formLowQualTag`.
- **MSI prediction for Ultima with the supplied model files:** null pointer; the model for Ultima uses
  AC/AG units, which the files do not have. `ms_table` is written, `jitter_params` is not.
- **Jitter fit with exactly two usable repeat counts in 7 to 15:** the regression throws "not enough
  data" and the run ends before `jitter_params` is written (seen with the classes, not in a full run).

### Dead or unused code
- `ReadCache.checkPosCheckAdjust` returns at once unless the pop distance already differs from the
  default, and nothing else changes it, so the "dynamic pop distance" never takes effect.
- `msi_min_map_quality`, `msi_max_sites_per_type` and `max_site_alt_contribution` are read and never
  used.
- The jitter analysis never increments its total read count, so the "20% of reads rejected" site filter
  cannot fire; the outlier pass in `buildJitterCountsTable` never marks an outlier.
- `-output_id` is not applied to the name of `bqr.tsv`.
- The comment in `PartitionReader` says Ultima does not mark duplicates; it does.

### BQR
- `mMaxIndex` is overwritten rather than maximised, so data at higher positions can be dropped at a
  purge or block end.
- A read that runs past the end of the 100K block is cut off.
- Short-fragment boundaries are not reset for unpaired reads and carry over from the last paired read.
- The mate end computed from the mate cigar ignores lengths before I, =, X and H.
- A probability of exactly 1 would print as `-0.00`.

### SBX
- `ConsensusReads.buildFromRead` shares the base and quality arrays of the read it copies, so
  finalising the consensus read also changes that duplicate.
- `SbxDuplicateCollapser.collapseKey` uses the lower orientation, which is forward for both strands of
  an unpaired read, so forward and reverse reads collapse together.
- `SbxRoutines.finaliseRead`, merging X into M: `remove(nextElement)` removes the first equal cigar
  element, which can be an earlier one.
- `SbxDuplexIndelBuilder.buildAdjustedSupplementaryCigar` drops the last element when there is no soft
  clip.
- `SbxConsensusTest.testMismatchRepeatMirroring` asserts nothing.

### Other
- htsjdk writes a record's attributes in file order until any attribute, the cigar, the bases, the
  qualities or the name is set, and sorted by tag afterwards. So in REDUX's output the attribute order
  of a read tells whether REDUX changed it. Not a defect, but it surprised us.
- `SAMRecord.equals` is reached through `List.remove` and `contains` in the SBX code; it decodes both
  records, which also changes how they are written.

## 4. What is faster and why

- Java calls samtools three times (sort, concatenate, index); the port sorts the small per-partition
  files itself, copies compressed blocks and indexes with htslib.
- Java spends about 44% of its time in BQR (one object with a map of lists per position) and about 14%
  decoding bases; the port counts into flat arrays and reads bases from the packed record.
- The index is built while the partitions are written (one index per partition file, merged with the
  offsets at which the files land in the BAM), so the finished BAM is not read again. On a 78 GB output
  on a network disk the re-read took 227 s. The index is not byte-identical to one from `samtools index`
  (chunks are merged differently); 3,307 region queries, `idxstats` and htsjdk reads gave the same
  results with both.
- The jitter analysis costs about one second where Java allocates for every read and builds every
  table twice.
- Memory is a few hundred MB against about 20 GB.

## 5. Timeline

- Day 2, later: first Linux build and runs (x86-64 cloud VM): whole genome in 26 minutes against
  Java's 50 at 16 threads, same comparison result as on macOS.
- Day 2: whole genome compared at 1 and 16 threads; the read cache's dependence on earlier partitions
  found and followed; index built while writing.
- Day 1: Java profiled; Illumina path, both unmapping passes, consensus, output; first identical BAM.
  BQR, jitter, SBX, Ultima, UMIs and the unit tests ported in parallel and merged.
