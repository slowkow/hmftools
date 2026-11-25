# REDUX in Rust: record of the port

A running record, kept for a later announcement and for reports to the Java maintainers. Status and
measurements are in [PLAN.md](PLAN.md).

## 1. How it was checked

Java REDUX 2.1 (this repository, `mvn -pl redux -am package`) and the port were run on the same input
with the same options. The BAMs were compared record by record with `redux-rs compare_bams`; the other
files with `cmp`.

Everything in this section is about the port in Java mode (`-java_compat`; until 2026-11-22 the only
mode). The default mode corrects the Java faults of section 3, so its output differs from Java's where
a fault had an effect: section 3.5 has the measured differences, section 6 the whole genome in both
modes.

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
| The same, after the second optimisation pass, local SSD | 16 threads | 1,780,126,609 | BAM and all four side files identical to Java's one-thread output; 16.9 min, 1.1 GB |
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

Unit tests: the Java tests of the ported classes are ported, and each correction has a test of its own
where one can be written without real data (211 Rust tests in all).

## 2. Where Java's output depends on the thread count

None of these is a difference between the port in Java mode and Java at one thread. The ids are those
of section 3. For each, "corrected mode" says whether the default mode gives one result at any thread
count.

- **Partitions (R15).** The genome is split into `2 x threads` partitions and duplicates are not
  matched across a boundary, so a few duplicate groups change with the thread count. The port splits
  the same way and, in Java mode, matches Java at the same thread count. *Corrected mode:* not
  removed; it needs a redesign of the partition stage (section 3.6). Its size is measured below.
- **BQR block grid and partition boundaries (R16).** The block grid starts at each partition's first
  read and a position near a boundary is filtered once per partition, so `bqr.tsv` changes with the
  thread count: on the small Illumina input Java's file at 64 threads differs from its file at 1
  thread in 854 of 2,816 rows. *Corrected mode:* one result. The window of R2 has no grid, and the
  positions next to a boundary are counted once with the reads of both partitions. With
  `-bqr_jitter_msi_only` the file is the same at 1, 4 and 64 threads. In a full run the
  recalibration sees the reads left after duplicate marking, so it inherits what R15 changes.
- **Unplaced reads (R17).** With several threads all workers append to one synchronised writer, so the
  order of the fully unmapped reads at the end of the BAM depends on timing. The port writes Java's
  one-thread order in both modes.
- **`duplicate_freq.tsv` labels (R18).** For group sizes above 10 the size written is the first one
  seen for its rounded bucket. Statistics are kept per thread and merged, so the label depends on
  which thread took which partition. In Java mode the port gives the one-thread label. *Corrected
  mode:* the label is the rounded size, which does not depend on any order.
- **`umi_coord_freq.tsv` (R19).** `UmiStatistics.merge` keeps the first thread's maximum for a row that
  exists in both, instead of the larger, so 14 of 2,038 rows differed between Java at 8 threads and at
  1. The port gives the one-thread result in both modes.
- **`ultima_stats.tsv` row order (R20)** is HashMap order and, with several threads, depends on the
  same. The port gives the one-thread order in both modes.
- **Duplicate groups in very deep regions (R21).** Each thread keeps one `ReadCache` for all the
  partitions it takes, and the cache's size-check counter (`mCheckSizeReadCount`, one check every
  10,000 reads) is not reset between partitions. When the cache holds 100,000 reads or more, the check
  decides when half of it is purged early, and reads that arrive later for a purged position form a
  new group. So the groups in a very deep region depend on how many reads the same thread cached
  before it: with one thread, on every earlier partition; with several threads, on which partitions
  that thread happened to take. Seen on a whole genome, on a 33 kb decoy contig holding 7.5 million
  reads: adding 31,000 unrelated reads on chr1 to the input changed the duplicate groups on the decoy
  at one thread. In Java mode the port gives Java's one-thread result at any thread count: a partition
  whose cache nears the threshold waits for the read counts of the earlier partitions. *Corrected
  mode:* the counter starts again with each region, so a region's groups depend on its own reads only
  (the contig gives the same records with and without the chr1 reads), and nothing waits. What is
  left is that a deep region cut by a partition boundary starts its count at the boundary; that is
  part of R15.
- **Interim BAM of unmapped reads (R22).** Written by several threads in timing order, then sorted by
  position; reads at the same position and strand keep the timing order. The port keeps the order of
  discovery in both modes.
- **UMI length (R37).** The length of a fixed-length UMI comes from the first read whose UMI is
  extracted, which with several threads is a matter of timing, in Java and in the port. It matters
  only when UMIs differ in length and `-umi_variable_length` is not given. Left.

The jitter analysis does not depend on the thread count in either mode.

Measured in the corrected mode, same input at different thread counts (worker threads limited to 6,
which does not change the partitions):

| Input, run | Partitions compared | Result |
|---|---|---|
| Illumina small, `-bqr_jitter_msi_only` | 1, 8 and 128 | all three files identical. With R16 switched off and R2 on, 1 row of `bqr.tsv` differs at 128 partitions; Java's own files at 1 and 128 partitions differ in 854 rows |
| Illumina slice, full run | 1, 16 and 128 | BAM and all side files identical |
| SBX tumor chr21, full run | 1, 16 and 128 | identical (at 128 partitions one boundary lies inside the reads) |
| Ultima tumor chr21 | 1 and 16 | identical |
| Deep contig | 1, 16 and 128 | identical |
| Illumina with duplex UMIs, full run | 1 and 16 | identical |
| The same | 1 and 128 | R15: 2 of 27.2 million records differ in the duplicate flag and 1 record is new (one group cut by a boundary); from them `bqr.tsv` differs in 120 of 4,608 rows, `duplicate_freq.tsv` in 2, `ms_table` in 1 and `umi_coord_freq.tsv` in 4 of 2,038 |

The whole genome was run in the corrected mode at 32 partitions only (section 6). In Java mode the
port's output at 32 partitions is identical to Java's at one, so on that input R15 changes nothing in
Java mode; the corrected mode at one partition was not run.

## 3. Java findings

Numbered R1 to R41. "Reproduced" means the port behaves the same way under `-java_compat`, so that its
output equals Java's; where Java stops with an exception, "reproduced as a failure" means the port
stops there too. "Not reproduced" means the port does not follow Java in either mode, because Java's
result there depends on thread timing or is a failure with nothing to compare. Line numbers are those
of the REDUX 2.1 source in this repository (`redux/src/main/java/com/hartwig/hmftools/redux/` unless
`hmf-common` is named).

Since 2026-11-22 the port corrects most of these by default. A finding marked **Corrected by default**
takes the corrected path unless `-java_compat` is given; what the corrected behaviour is and what it
changes is in section 3.5, and what was left and why in section 3.6. Every statement of sections 1, 2
and 4 about equality with Java is about the port run with `-java_compat`.

### 3.1 Wrong result

- **R1. `bqr/BqrRegionReader.java:524`** `mMaxIndex = posIndex` is set whenever a position is first
  touched, also when the index is lower than the highest so far. Reads reach the recalibration
  slightly out of order (they come from the read cache, before the sorted writer), and a read with a
  deletion or skip makes high positions before a later read makes lower ones. When the block ends,
  `buildQualityCounts` (line 222) stops at `mMaxIndex` and the positions above it are dropped.
  Reproduced. **Corrected by default.**
- **R2. `bqr/BqrRegionReader.java:154-155, 406`** The positions are held in blocks of 100,000 bases
  that start at the first read after the previous block. A read that runs past the end of its block is
  cut off (`afterEnd`, line 406), and those positions get reads only from the next block on. A read
  that arrives out of order with a start before the current block ends the block and starts one at the
  lower position (line 137-160), so the positions of the abandoned block are counted in two parts,
  each with its own alt and indel filters. Reproduced. **Corrected by default.**
- **R3. `bqr/BqrRegionReader.java:298-299`** `setShortFragmentBoundaries` returns for an unpaired read
  before it resets the two boundaries, so an unpaired read is clipped by the boundaries of the last
  paired read. It needs an input that mixes paired and unpaired reads. Reproduced. **Corrected by
  default.**
- **R4. `hmf-common` `bam/CigarUtils.java:201`** `getReadBoundaryPosition` over a cigar string counts
  the elements D, M, S and N. The lengths of `=` and `X` elements are dropped, so for a mate cigar
  written with `=` and `X` the mate's end is its start minus one. It is used for the short-fragment
  boundary of the recalibration (`bqr/BqrRegionReader.java:324`), for the fragment coordinate of a
  reverse mate (`duplicate/FragmentCoords.java:183`), and for the coordinates taken from a
  supplementary alignment's cigar; `calcCigarAlignedLength` has the same list. The lengths before I
  and H are dropped too, which is right. No effect on bwa output, which uses M. Reproduced.
  **Corrected by default.**
- **R5. `bqr/BqrRegionReader.java:115-116, 214`** `onRegionComplete` leaves `mCurrentRefSequence` in
  place, and `initialise` returns at once for a region that is not on a human chromosome, leaving the
  reader active. The reference block of one region is then the current block of the next, and a read
  of the next region whose start falls inside its positions is counted against the bases of the wrong
  chromosome. With reads in genome order and the standard regions the positions do not coincide; with
  `-specific_regions`, or when a thread takes a partition that starts inside a chromosome, they can.
  Found while writing the corrections. Reproduced. **Corrected by default.**
- **R6. `consensus/ConsensusReads.java:252-253`** `buildFromRead` gives the new record the base and
  quality arrays of the read it copies. `SbxRoutines.finaliseRead` then changes the consensus read's
  bases and qualities and with them those of that duplicate, which is written with them when it is
  written from its decoded fields. Reproduced. **Corrected by default.**
- **R7. `duplicate/SbxDuplicateCollapser.java:460`** `collapseKey` takes the orientation of the lower
  position. For an unpaired read the lower position of a reverse read is its 3' end, with the
  opposite orientation (`duplicate/FragmentCoords.java:218`), so the key is "forward" on both strands
  and forward and reverse reads within the collapse distance are put in one duplicate group. The class
  comment (line 42) and the fragment key (`FragmentCoords.java:84`) keep the strands apart; the last
  block of `UnpairedReadDuplicateTest.testSbxDuplicateGroupCollapsing` passes only when they are
  collapsed. In a run the two strands of one span leave the read cache at different times, so they
  meet in the collapser only at the end of a region or when the reads are shorter than the pop
  distance. Reproduced. **Corrected by default; intent uncertain** (3.5).
- **R8. `consensus/SbxRoutines.java:822`** Merging `X` elements into `M`: `cigarElements.remove(nextElement)`
  removes the first element equal to the next one, which is an earlier element when one has the same
  length and operator (`5M2D3M1X5M` becomes `2D9M5M`). Reproduced. **Corrected by default.**
- **R9. `consensus/SbxDuplexIndelBuilder.java:356`** `buildAdjustedSupplementaryCigar` assumes a soft
  clip. For a read with supplementary data and no soft clip it drops the read's last cigar element and
  appends what is left of the supplementary cigar. Reproduced. **Corrected by default.**
- **R10. `duplicate/UmiGroupBuilder.java:485`** `splitUmi` takes the second half of a duplex UMI from
  one character after the start of the delimiter, whatever the delimiter's length, so with a
  delimiter of two characters the second half begins with the delimiter's second character.
  Reproduced. **Corrected by default.**
- **R11. `jitter/MicrosatelliteSiteData.java:158`** `mTotalReadCount` is never incremented
  (`addReadToStats`, line 100), so the rejected fraction is `1 - passing / 0.0`, minus infinity, and
  the README's rule "20% or more considered reads are rejected" never excludes a site. Reproduced.
  **Corrected by default.**
- **R12. `jitter/MsJitterAnalyser.java:144-193`** `buildJitterCountsTable` builds each table in a loop
  "until no outlier is found", but nothing sets `outlierFound` or adds to `outliers`, so the second
  pass repeats the first. Reproduced as a single pass with the same result. Left (3.6).
- **R13. `jitter/MsJitterConfig.java:81-83`** `msi_min_map_quality`, `msi_max_sites_per_type` and
  `max_site_alt_contribution` are read and never used. The reads of the jitter analysis are filtered
  by the minimum of the recalibration (`write/BamWriter.java:141`). Reproduced. **The first is
  corrected by default**; the other two are left (3.6).
- **R14. `duplicate/ReadCache.java:539`** `checkPosCheckAdjust` returns whenever the pop distance is
  the default, and nothing else changes it, so the distance is never shortened when the cache holds
  100,000 reads or more ("creep in 10% of the usual pop distance check", line 550). Reproduced.
  **Corrected by default; intent uncertain** (3.5).
- **R18. `common/Statistics.java:80`** A row of `duplicate_freq.tsv` for group sizes above 10 stands
  for a range of sizes (`roundFrequency`), but it is labelled with the first size seen for it
  (`new DuplicateFrequency(duplicateCount, ...)`), not the rounded size. The label then depends on the
  order of the groups, and with several threads on which thread took which partition. Reproduced at
  one thread; with more threads the port gives the one-thread label. **Corrected by default.**

### 3.2 Nondeterminism: output that depends on the thread count or the run

Section 2 says for each of these what the corrected mode gives at different thread counts.

- **R15. Partitions (`PartitionThread.splitRegionsIntoPartitions`, `ReduxApplication.java:224`).** The
  genome is split into `2 x threads` partitions and a read belongs to the partition that holds its
  alignment start. Reads of one duplicate group whose starts lie on both sides of a boundary are not
  matched. Reproduced: the port splits the same way. Left (3.6).
- **R16. `bqr/BqrRegionReader.java:154`** The block grid of the recalibration starts at the first read
  of each partition, and a position near a partition boundary is counted once with the reads of each
  partition, each time through the alt and indel filters. So `bqr.tsv` changes with the thread count.
  Reproduced. **Corrected by default**, together with R2.
- **R17. `write/BamWriterSync.java:40`** With several threads all workers append the fully unmapped
  reads to one synchronised writer, so their order at the end of the BAM depends on timing. Not
  reproduced: the port writes Java's one-thread order at any thread count.
- **R19. `duplicate/UmiStatistics.java:69-74`** `merge` takes the maximum group of a row from the
  other thread only when its own row is new, and keeps its own otherwise, instead of the larger of the
  two. Which partitions a thread takes depends on timing. Seen: 14 of 2,038 rows of
  `umi_coord_freq.tsv` differed between Java at 8 threads and at 1. Not reproduced: the port gives the
  one-thread result.
- **R20. `consensus/UltimaStats.java:23`** The rows of `ultima_stats.tsv` are written in HashMap order,
  which with several threads depends on the order in which the threads' maps were merged. The
  one-thread order is reproduced; at more threads the port gives the one-thread order.
- **R21. `duplicate/ReadCache.java:511, 570-577`** Each thread keeps one `ReadCache` for all the
  partitions it takes. The size-check counter (`mCheckSizeReadCount`, one check every 10,000 reads)
  and the last counted size are reset only by `clear()`, which nothing calls; `evictAll` at the end of
  a region leaves them. When the cache holds 100,000 reads or more the check decides when half of it
  is purged early, and reads that arrive later for a purged position form a new group. So the groups
  of a very deep region depend on how many reads the same thread cached before it: at one thread on
  every earlier region, at several threads on which partitions that thread happened to take. Seen on
  a whole genome on a 33 kb decoy contig holding 7.5 million reads: adding 31,000 unrelated reads on
  chr1 changed the duplicate groups on the decoy at one thread. Reproduced as the one-thread result at
  any thread count. **Corrected by default.**
- **R22. `write/FileWriterCache.java:92`** The interim BAM of unmapped reads is written by several
  threads in timing order and then sorted by position; reads at the same position and strand keep the
  timing order. Not reproduced: the port keeps the order of discovery.
- **R37. `duplicate/UmiConfig.java:151-157`** The length of a fixed-length UMI is taken from the first
  read name whose UMI is extracted. With several threads that is whichever partition gets there
  first. It matters only when the UMIs differ in length and `-umi_variable_length` is not given. The
  port has the same dependence in both modes. Left (3.6).

- **R41. Results depend on the order of the input BAMs** (`BamReader.java`, `consensus/IlluminaRoutines.java:230-255`,
  2026-11-25). Several input BAMs (`-input_bam a,b,c`, one per lane) are merged by alignment start, and
  reads that start at the same position come in the order the files were given. That order reaches the
  results where the code keeps the first of equals: the consensus base of a position at which two bases
  have the same quality total and neither is the reference (always so in soft-clipped bases, which have
  none) is that of the first read, and so are other choices between equals. Observed on Google Batch: a
  35x tumor from three lanes, the same lane BAMs given in two orders (oncoanalyser listed them in the
  order the alignment tasks finished), 16 threads, corrected mode: 7,751 of 1.59 billion records
  differed (60 in the duplicate flag, 7,441 in bases and 273 in qualities: consensus reads, the only
  records whose bases REDUX rewrites), 16 and 24 records only in one; `duplicate_freq.tsv` differed in 8
  of its rows; the 4-lane normal, listed in the same order both times, was byte-identical. On the
  two-lane test input of section 1 the two orders give the same records (Java mode) except the order of
  32 unplaced reads. Reproduced. Not yet corrected; the slowkow/oncoanalyser fork now lists the lane
  BAMs by file name, so a pipeline run gives one result.

### 3.3 Crash

- **R23. `consensus/BaseBuilder.java:148`** Ultima with `-form_consensus`: `determineBaseAndQual`
  returns an invalid base for Ultima ("TODO: decide for Ultima") wherever the reads of a duplicate
  group disagree, and the BAM encoder then throws "Bad base passed to charToCompressedBaseLow". Two
  reads at one position that differ in one base reproduce it. With `-no_write_bam` Java completes,
  having passed the invalid base to the recalibration. The port stops at the same read with a message
  in both modes, also with `-no_write_bam`. Left (3.6).
- **R24. `duplicate/UmiGroupBuilder.java:561`** `-umi_base_diff_stats`: `captureUmiGroupStats` casts a
  `ReadInfo` to `SAMRecord`, so the first fragment without a duplicate ends the run with a
  ClassCastException. The option cannot complete on ordinary data. Reproduced as a failure.
  **Corrected by default.**
- **R25. `duplicate/UmiGroupBuilder.java:484-485`** A duplex UMI without the delimiter: `indexOf`
  gives -1 and `substring` throws. Reproduced as a failure. **Corrected by default.**
- **R26. `hmf-common` `sequencing/UltimaBamUtils.java:53`** Ultima with `-ultima_drop_qual_tags` and
  the recalibration on: the `t0` attribute is dropped before the recalibration reads it, a
  NullPointerException in `extractT0Values`. Reproduced as a failure. **Corrected by default.**
- **R27. `consensus/UltimaRoutines.java:106-131`** An Ultima read without `t0`, or with one shorter
  than the read: NullPointerException or index error in `formLowQualTag`, and the same in the
  recalibration for a BAM that was written without the attribute. Reproduced as a failure with a
  message; the recalibration path of the port ended with a panic. **Corrected by default** (one
  message, status 1).
- **R28. `ms_model/MsModelCalculator.java:350`** MSI prediction for Ultima with the supplied model
  files: the model for Ultima uses the AC and AG units, which the files do not hold;
  `mRepeatUnitCountErrorRates.get` is unboxed without a null test and the run ends before
  `jitter_params` is written (`ms_table` is written). The next line tests the coefficients for null.
  Reproduced as a failure. **Corrected by default.**
- **R29. `jitter/JitterModelFitter.java:170-172`** With exactly two repeat counts in 7 to 15 that have
  a row, both with 20,000 reads or more, the test `numRepeatCountsOverMin > 1` passes and the
  regression throws "not enough data" (Commons Math asks for more rows than two predictors), and the
  run ends before `jitter_params` is written. Seen with the classes, not in a full run. Reproduced as a
  failure. **Corrected by default.**
- **R38. `duplicate/UmiStatistics.java:287-296`** `-umi_base_diff_stats` on an input in which no UMI
  is recorded: `writeUmiBaseFrequencyStats` walks a table that was never made. Not reproduced: the
  port writes the header alone, in both modes.
- **R39. `duplicate/UmiConfig.java:163`** A read name shorter than the UMI length: `substring` throws.
  Not reproduced: the port takes the whole name, in both modes.
- **R40. `hmf-common` `bam/SamRecordUtils.java:134`, `bam/CigarHandler.java:82`** A `UT` attribute
  other than the three consensus type names, or a `P` cigar element, ends the run with an exception.
  In the jitter analysis the port reads such a type as NONE and skips the element, in both modes; in
  the recalibration it stops as Java does.

### 3.4 Oddity

- **R30. `bqr/BaseQualRecalibration.java:258`, `common/Statistics.java:139`,
  `jitter/MsJitterAnalyser.java:103, 135`** `-output_id` is applied to the BAM and to the UMI and
  Ultima statistics files (`ReduxConfig.formFilename`), not to `bqr.tsv`, `duplicate_freq.tsv`,
  `jitter_params.tsv` and `ms_table.tsv.gz`, whose names come from `hmf-common`. Reproduced. Left
  (3.6).
- **R31. `bqr/BaseQualRecalibration.java` (`calcRecalibratedQual`)** A probability of exactly 1 gives
  `-10 * log10(1) = -0.0`, printed as `-0.00`. Not seen in any output. Reproduced. **Corrected by
  default.**
- **R32. `PartitionReader.java:100`** The comment says Ultima does not mark duplicates; it does.
- **R33. `SbxFinaliseReadTest.testMismatchRepeatMirroring`** (test, line 196) asserts nothing.
- **R34. htsjdk** writes a record's attributes in file order until any attribute, the cigar, the
  bases, the qualities or the name is set, and sorted by tag afterwards. So in REDUX's output the
  attribute order of a read tells whether REDUX changed it. Not a defect, but it surprised us.
  Reproduced in both modes, so that the BAMs of the two modes can be compared record by record.
- **R35. `SAMRecord.equals`** is reached through `List.remove` and `contains` in the SBX code; it
  decodes both records, which also changes how they are written. Reproduced in both modes, for the
  same reason.
- **R36. `bqr/BaseQualityResults.java:38-39`** `addBaseQualityRegionCounter` adds the reader's running
  totals of reads used and alts filtered each time a block completes, so the two numbers of the debug
  line count the early blocks many times. Reproduced. **Corrected by default** (log only).

### 3.5 Corrections

24 of the findings above are corrected by default. Each corrected place in the code asks
`compat::java("<id>")` (`src/compat.rs`) and keeps the Java path beside the corrected one.

- `-java_compat`, or `REDUX_JAVA_COMPAT=1` in the environment: Java behaviour for every finding.
- `-java_compat_except R7,R2` (or `REDUX_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
  those listed. The effects below were measured this way, one finding at a time, against the stored
  Java output of the same run (to which the port in Java mode is identical).

Inputs of the measurements: "Illumina" is the realigned slice of 25.8 million reads (27.7 million
records out); "SBX" and "Ultima" are the chr21 slices of GIAB HG008 tumor (23.8 and 21.2 million
records); "deep contig" is a slice of one decoy contig holding 7.5 million reads in 33 kb, alone and
together with 31,000 reads on chr1; "UMI" is the Illumina slice with simulated UMIs. BAM effects are
counts of records by `redux-rs diff_bams`: changed in a field, gone, or new. Rows of `bqr.tsv` are
counted without the header. "Intent uncertain" marks a correction where the Java code allows more
than one reading; the reading taken and the alternatives are given. These should be confirmed with
the REDUX authors.

**Base quality recalibration**

| Id | Corrected behaviour | Effect |
|---|---|---|
| R1 | the max index is the highest position made | `bqr.tsv`: Illumina 440 of 2,816 rows, 37,241 more bases counted of 3.24 billion; SBX (a 10 Mb region) 109 of 2,417 rows; Ultima 645 of 3,072 rows |
| R2 | the positions are held in a window that moves with the reads: no read is cut off, a late read does not restart anything, and every position is filtered once with all its reads. The window starts 2,500 bases before a read, the distance at which Java purges; a base that arrives for a position already purged is not counted (Java makes the position again and never counts it) | `bqr.tsv`: Illumina 1,357 of 2,816 rows, 739,355 more bases (0.02%); SBX region 819 of 2,417 rows, 243,423 more bases (0.03%); Ultima 814 of 3,072 rows. With the window the counts no longer depend on the order in which the reads arrive |
| R16 | with R2: the positions next to a partition boundary are kept back by both partitions and counted once with the reads of both; a partition keeps back its first 100,000 positions for the reads of the partition before | section 2: `bqr.tsv` is the same at 1, 4 and 64 threads with `-bqr_jitter_msi_only` |
| R3 | an unpaired read has no short-fragment boundaries | none on any input (none mixes paired and unpaired reads) |
| R4 | `=` and `X` elements of a cigar string count as `M` does, in the mate's end, the fragment coordinates and the aligned length | none on any input (bwa writes M) |
| R5 | the reference block ends with its region, and a region off the human chromosomes leaves the reader inactive | none on any input |
| R31 | a probability of 1 prints as `0.00` | none on any input |
| R36 | the two totals of the debug line are added once | debug log only |

**Duplicates and consensus**

| Id | Corrected behaviour | Effect |
|---|---|---|
| R6 | a consensus read copied from one read has its own bases and qualities | SBX: 3,410 records changed (qualities 3,410, bases 722), none gone or new. Illumina: not reached |
| R7 | **intent uncertain.** The collapse key of an unpaired read takes the read's own orientation, as the fragment key does, so forward and reverse reads are collapsed separately (the class comment). Alternative: Java's behaviour is meant, as the last block of its unit test expects; the strands would then be collapsed only when they happen to leave the cache together | SBX: 1,201 records changed (1,000 duplicate flags), 399 gone, 297 new, 102 fewer records; 14 of 30 rows of `duplicate_freq.tsv`. On a 10 Mb region of the same input: no record |
| R8 | the element that was merged is the one removed | none on the SBX inputs |
| R9 | a read with supplementary data and no soft clip keeps its cigar | none on the SBX inputs |
| R14 | **intent uncertain.** `checkPosCheckAdjust` returns early only when the distance is the default and the cache is under the threshold, so the pop distance goes down by 10% of the default for each 100,000 cached reads and back up below the threshold, as the rest of the method does. Alternative: the early return is a deliberate switch-off and the early purge alone is meant to bound the cache | deep contig: 6,147 records changed (2,688 duplicate flags), 830 gone, 654 new, 176 fewer records (fewer groups split by an early purge). Illumina, SBX: none |
| R21 | the size-check counter and the last counted size start again with each region, as after `clear()`; no partition waits for the read counts of earlier ones | deep contig alone: none. With the reads on chr1: 400 records changed (171 duplicate flags), 68 gone, 77 new, and the records on the contig are then those of the contig alone (in Java 171 flags differ between the two inputs). Illumina, SBX: none |
| R18 | the label of a row is the rounded size | Illumina 1 of 10 rows; SBX (10 Mb region) 19 of 28 rows; labels only, the counts are the same |

**UMIs**

| Id | Corrected behaviour | Effect |
|---|---|---|
| R10 | the second half of a duplex UMI starts after the whole delimiter | no input has such a delimiter; unit test |
| R24 | a fragment without duplicates counts as a group of its one read (`new DuplicateGroup(umi, read, null)`, which the Java line intends) | `-umi_base_diff_stats` on the UMI input completes where Java stops; BAM identical to the run without the option, `umi_coord_freq.tsv` differs from it in 4 of 2,042 rows; `umi_edit_distance.tsv` (11 rows) and `umi_nucleotide_freq.tsv` (8 rows) are written |
| R25 | a UMI without the delimiter is no duplex match; one warning names the delimiter | `-umi_duplex` on the single-UMI input completes where Java stops |

**Jitter and MSI**

| Id | Corrected behaviour | Effect |
|---|---|---|
| R11 | the total is every read of the site, passing or rejected, so a site with 20% or more rejected reads is left out | Illumina: 864 of 2.38 million site rows become `realVariant`, 51 of 89 rows of `ms_table` and 4 of 7 rows of `jitter_params` change. SBX: 9,996 of 4.75 million site rows, 104 of 106 and 7 of 14. Ultima: 124 of 174 and 6 of 14 |
| R13 | `-msi_min_map_quality` is the lowest mapping quality of a read of the jitter analysis (default 50, as before) | none at the default |
| R28 | a unit and count without an error rate is left out of the prediction, as one without coefficients is; one warning | the Ultima run with the supplied model files completes: `msi_prediction.tsv` holds -1 (no unit of the Ultima model is in the files) and `jitter_params.tsv` is written |
| R29 | two rows give the line through them; it is used when it rises, as with more rows | none on any input (no table has exactly two such rows) |

**Ultima**

| Id | Corrected behaviour | Effect |
|---|---|---|
| R26 | with `-ultima_drop_qual_tags` the `t0` and `tp` attributes stay on a read until the recalibration has counted it and are taken off as it is written | the run completes where Java stops. Its BAM is identical to Java's with `-bqr_disable` added, and its `bqr.tsv` to Java's without the option |
| R27 | one error line naming the read and status 1, also in the recalibration | message |

All corrections together, default against `-java_compat`:

| Input | BAM | Side files |
|---|---|---|
| Illumina, one BAM and two lane BAMs | identical | `bqr.tsv` 1,365 of 2,816 rows; `duplicate_freq.tsv` 1 label; `jitter_params` 4 of 7 rows; `ms_table` 52 of 89 rows |
| Illumina with duplex UMIs | identical | `bqr.tsv` 1,945 of 4,608 rows; `duplicate_freq.tsv` 1 label; `jitter_params` 6 of 21; `ms_table` 147 of 267; UMI statistics identical |
| SBX tumor | 4,611 of 23.8 million records changed (qualities 4,577, bases 1,272, duplicate flag 1,000), 399 gone, 297 new; unplaced reads the same | `bqr.tsv` 1,060 of 2,493 rows; `duplicate_freq.tsv` 30 of 30; `jitter_params` 7 of 14; `ms_table` 104 of 106 |
| SBX normal | 1,003 of 11.0 million records changed (qualities 999, duplicate flag 226), 96 gone, 16 new | `bqr.tsv` 1,028 of 2,406 rows; `duplicate_freq.tsv` 25 of 29; `jitter_params` 7 of 14; `ms_table` 124 of 124 |
| Ultima tumor | identical | `bqr.tsv` 814 of 3,072 rows; `jitter_params` 6 of 14; `ms_table` 124 of 174; `ultima_stats.tsv` identical |
| Deep contig | 6,147 of 7.56 million records changed, 830 gone, 654 new | `duplicate_freq.tsv` labels |

Unit tests: every correction that can be shown without real data has a test of the corrected path
beside a check of the Java path (R1 to R5, R7 to R11, R13, R14, R16, R18, R21, R24, R25, R26, R28, R29,
R31). The Java tests that pass only with a fault run inside `compat::with_java(true, ...)`: the last
block of `testSbxDuplicateGroupCollapsing` (R7) and the ported checks of the failures. R6 and R27 have
no unit test: R6 needs a duplicate group through the partition reader and is covered by the SBX runs,
R27 ends the process.

### 3.6 Left as in Java, and why

- **R12, outlier pass.** Nothing in the code says what an outlier is: the loop holds an empty set and a
  flag that is never set, and `max_site_alt_contribution` (default 0.8, "Max percentage a single alt
  site can contribute") is not used anywhere. The README does not describe an outlier rule. A rule
  could be guessed from the option's description (leave out a site that holds more than 80% of the
  reads of one jitter in a table row), but that would be inventing behaviour.
- **R13, `msi_max_sites_per_type` and `max_site_alt_contribution`.** The first belongs to choosing
  sites, which REDUX no longer does (the site file is made by the separate site finder); the second
  to the outlier pass above. Both are still read and checked.
- **R15, duplicates across partition boundaries.** A read belongs to the partition of its alignment
  start, but a duplicate group is defined by fragment coordinates: the unclipped 5' position of a
  forward read lies before its start, that of a reverse read after it, and UMI jitter merging and the
  SBX collapse look a few bases further. To give one result at any thread count each partition would
  have to read the reads within that distance on both sides of its boundaries, let them take part in
  grouping without writing them, and write only the records (and the consensus read) that start
  inside it, with the statistics counted by one of the two partitions. That changes how the read
  cache pops, the bounds of the sorted writer, the UMI and SBX collapsing and the statistics, and the
  early purge of a deep region (R21) would have to depend on the cache content alone, not on a count
  of reads since the region began. It is a redesign of the partition stage and was not done. Measured
  size of the dependence in the corrected mode: section 2.
- **R23, Ultima consensus.** Java has no rule for a base on which the reads of a group disagree, and
  nothing near it offers one: a consensus read for Ultima would also need consensus values for the
  `t0` and `tp` attributes that the recalibration and the low-quality attribute are made from. The
  port stops with one message that names the read and tells the user to run Ultima without
  `-form_consensus`. It also stops with `-no_write_bam`, where Java completes with an invalid base.
- **R30, `-output_id`.** The four files without the id are found by name by the tools that read them
  (SAGE reads `bqr.tsv` and `jitter_params.tsv` as `<sample>.redux.<type>.tsv`), and the four are
  treated alike, so this looks like a choice, not an oversight. Applying the id would break the
  hand-over.
- **R37, UMI length.** One result at any thread count would need the UMI length before the partitions
  start (from the first read of the input, say), which is not what one thread of Java does either
  (it takes the first read whose UMI is extracted). It matters only for inputs that should be run
  with `-umi_variable_length`.
- **R17, R19, R20, R22, R38, R39, R40** are not reproduced in either mode (see each finding): the port
  already gives one result where Java depends on timing, or completes where Java fails.
- **R32, R33** are a comment and a test. **R34, R35** are properties of htsjdk that the port follows
  in both modes so that its records can be compared with Java's.

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
- **Port bug, found on the first Google Batch run (2026-11-23) and fixed:** until then, a record that
  followed a block ended early (because the record did not fit) started, in the index, at the end of
  the previous block rather than at the start of the next one. htslib reads such offsets; htsjdk
  refuses to seek there (`Invalid file pointer`), so Java BamTools failed on the REDUX BAM, and any
  Java tool querying it by region could. On a 2M-read-pair BAM, 1,564 of 66,021 chunk starts were
  affected; the earlier htsjdk checks had not queried any of them. The writer now does what htslib's
  own indexed writer does (`bgzf_flush_try`, then `hts_idx_amend_last`). The BAM is byte-identical; the
  index now gives htsjdk's answers for every 100 kb window, the same records as `samtools index`
  (35,264 queries, 28.3 million records), and the unit test
  `out::tests::index_offsets_never_point_at_a_block_end` checks it (28 bad offsets before the fix).
- The jitter analysis costs about one second where Java allocates for every read and builds every
  table twice.
- Memory is a few hundred MB against about 20 GB.
- A later pass over the profile cut the port's own CPU by about a third on a one-thread run (vendor
  chr21 slice: 98.9 s to 63.6 s of user time, 1,313 to 765 billion instructions), with identical output:
  - the fragment key of each read is built in one allocation, and the read cache keeps each key once
    in a table hashed with Java's hash, instead of copying and re-hashing it (18% of instructions);
  - a faster allocator, mimalloc, is the default build (13%);
  - the recalibration loop works out the positions to count once per aligned block rather than testing
    each base (12%);
  - reference bases: the FASTA is read at offsets worked out from the index rather than a byte at a
    time through htslib, the previous block is kept, and ranges are copied by block (about 5%);
  - attributes are found by walking the record's attribute block, not through htslib's lookup (2%).
- A third pass took another 6% (765 to 720 billion instructions on the same run), identical output:
  a consensus path for the usual case of reads that all cover every base, bases decoded two at a
  time (3%); the sorted writer holds a list of records per start position rather than one tree of
  whole records (2%); string attributes are skipped with a vectorised search (1%).
- The partition files are now appended to the output BAM while later partitions are still running,
  and each is deleted once appended, so less temporary disk is needed. On the whole genome (16
  threads, external USB SSD) this gave no time saving (17.2 against 16.9 minutes): that disk writes
  about 300 MB/s of compressed data, and the second copy of the output is bound by it whenever it
  happens. A disk with headroom should gain up to half of the joining time; not measured.
- What remains is mostly output compression, about half of the CPU. It is at the fastest setting that
  keeps Java's file size; `-compress_level 3` saves about 9% of CPU for a 4.8% larger BAM and
  `-compress_level 1` about 21% for 7% larger (measured before the later changes above).

## 5. Timeline

- Day 5 (2026-11-22): the Java faults numbered (section 3) and corrected by default behind
  `-java_compat`; each correction measured on its own; whole genome in both modes.
- Day 4: third optimisation pass (consensus, sorted writer, attributes); partition files joined while
  the partitions run.
- Day 3: second optimisation pass (allocation, recalibration loop, reference access); a regression
  script reruns all 88 stored Java baselines.
- Day 2, later: first Linux build and runs (x86-64 cloud VM): whole genome in 26 minutes against
  Java's 50 at 16 threads, same comparison result as on macOS.
- Day 2: whole genome compared at 1 and 16 threads; the read cache's dependence on earlier partitions
  found and followed; index built while writing.
- Day 1: Java profiled; Illumina path, both unmapping passes, consensus, output; first identical BAM.
  BQR, jitter, SBX, Ultima, UMIs and the unit tests ported in parallel and merged.

## 6. Whole genome in both modes

The whole-genome input of the stored Java baselines (vendor alignments, 86 GB, 1.56 billion reads),
`-form_consensus`, the 32 partitions of `-threads 16`, worker threads limited to 6 (`REDUX_WORKERS=6`;
the machine was shared), external SSD, 2026-11-22.

Each mode was run twice; the second pair ran one after the other with the resident memory sampled
every 5 seconds. CPU is user plus system time.

| Run | Mode | Wall | CPU | Peak memory |
|---|---|---|---|---|
| first pair | `-java_compat` | 2,202 s | 5,704 s | 1,147 MB |
| first pair | default | 1,666 s | 5,647 s | 1,280 MB |
| second pair | `-java_compat` | 1,729 s | 5,433 s | 1,123 MB |
| second pair | default | 1,537 s | 5,373 s | 1,159 MB |

Java mode, first run: the BAM is identical to Java's one-thread BAM (1,780,126,609 records), and
`bqr.tsv`, `duplicate_freq.tsv`, `jitter_params.tsv` and `ms_table.tsv.gz` are identical to Java's.

Default against `-java_compat` (`redux-rs diff_bams`):

- **BAM.** 1,780,126,056 records against 1,780,126,609. The first 1,706,024,006 records are
  identical. After that, on the one very deep decoy contig, 6,867 records are changed (duplicate flag
  2,961, qualities 3,786, bases 2,468, attributes 3,904), 1,171 are gone and 618 are new: the effect of
  R14 and R21, as on the slice of that contig in section 3.5. Nothing else in the genome changes. The
  66,596,894 unplaced records are the same.
- **`bqr.tsv`:** 1,199 of 2,304 rows (R1, R2). **`duplicate_freq.tsv`:** 19 of 19 rows (the labels,
  R18, and the counts of the deep contig). **`jitter_params.tsv`:** 6 of 7 rows, and 121 rows of
  `ms_table.tsv.gz` (R11).

**Time.** The default mode took 24% and 11% less wall time for 1% less CPU time. In Java mode a
partition whose read cache nears 100,000 reads waits for all earlier partitions, to reproduce the
state of Java's one read cache (R21); in the default mode nothing waits. The wall times of the first
pair were taken while other jobs ran on the machine.

**Memory.** The peak of the default mode was 12% and 3% (134 and 36 MB) above that of Java mode. The
corrections do not use more memory for the same work: the two memory profiles are the same until the
last minute (about 0.9 GB), both peak in a spike of about 230 MB while the very deep contig is
processed, and that contig alone at one thread peaks at 625 and 629 MB in the default mode against
648 and 649 MB in Java mode (R14 shortens the pop distance). The difference is in what else runs
during the spike: in Java mode the partition with the deep contig has waited until every earlier
partition was complete, in the default mode it runs beside the other partitions. Bringing the peak
down to Java mode's would mean making the deep partition wait again, which is what costs Java mode its
wall time; this was not done. The peak also varies between runs of one mode by as much as the modes
differ in the second pair.
