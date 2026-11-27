# BamTools metrics in Rust: record of the port

A record of how `bam-tools-rs` was checked against Java BamTools 1.6.2 `BamMetrics` (this repository,
`mvn -pl bam-tools -am package`), the faults found in the Java code, and timings. Options and outputs
are in [README.md](README.md).

## 1. How it was checked

All comparisons run the port with `-java_compat`, 2 threads unless stated, against Java 1.6.2 (`-Xmx16g`,
2 threads), with the oncoanalyser command line (`-sample -bam_file -ref_genome -ref_genome_version 38
-driver_gene_panel -ensembl_data_dir -threads -output_dir`) plus the options named. A regression script
reruns every Rust case of the table below against the stored Java outputs (19 checks; all pass with the
final binary).

How each output is compared:

| File | Comparison | Why |
|---|---|---|
| summary, coverage, frag_length, flag_counts, gene_coverage, exon_coverage, off_target_overlaps | byte-identical | deterministic in Java |
| partition_stats | identical apart from the ProcessTime column, rows sorted | ProcessTime is a timing; Java writes rows in the order its threads finish partitions (BT9); the port writes partition order |
| target_regions, off_target_enriched | identical with rows sorted | written in thread order by Java (BT9) |

Results (Java mode); "identical" means every output file as in the table above:

| Run | Input | Options | Result |
|---|---|---|---|
| small_t | a small private tumor BAM (chr21, chr22, parts of chr1-10; 60.3M reads) | oncoanalyser command | 7 files identical |
| small_n | small normal | oncoanalyser command | 7 files identical |
| sbx_t | HG008 SBX tumor, chr21, after REDUX (18.6M reads; consensus reads, 2.86M dual-strand; 2.38M unplaced unmapped) | oncoanalyser command | 7 files identical |
| sbx_n | HG008 SBX normal, chr21, after REDUX | oncoanalyser command | 7 files identical |
| sbx_n_cram | the same as CRAM | oncoanalyser command | 7 files identical (and equal to the BAM's) |
| ultima_n | HG008 Ultima normal, chr21, after REDUX (37M reads, 2.56M unplaced unmapped) | oncoanalyser command | 7 files identical |
| ultima_t | HG008 Ultima tumor, chr21, after REDUX | oncoanalyser command | 7 files identical |
| v_quals | small tumor | `-specific_chr chr21 -map_qual_threshold 30 -base_qual_threshold 25 -max_coverage 100` | identical |
| v_nodup_zero | small tumor | `-specific_chr chr22 -no_duplicates -exclude_zero_coverage` | identical |
| v_regions | small tumor | `-specific_regions "chr21:20000000-30000000;chr22:30000500-31500000" -partition_size 300000` | identical |
| v_target_only | small tumor | `-regions_file testdata/synthetic_targets.bed -only_target -write_off_target -off_target_frag_overlap_threshold 5` | 9 files identical (target_regions, off_target_overlaps, off_target_enriched too) |
| v_target | small tumor | the same BED, `-specific_chr "chr21;chr22" -write_off_target -off_target_frag_overlap_threshold 3 -max_coverage 300` | 10 files identical |
| v_threads1 | small tumor | `-threads 1 -specific_chr chr21 -output_id run1 -partition_size 2500000` | identical; partition_stats rows in the same order without sorting |
| e_partition | small tumor | `-specific_chr chr22 -partition_size 100000` | identical |
| e_tsv | small tumor | `-regions_file testdata/synthetic_targets.tsv` (TSV with a header, a nested region) `-only_target -write_off_target -off_target_frag_overlap_threshold 2` | identical |
| e_noversion | small tumor | no `-ref_genome_version`, `-specific_chr chr21` (BT10) | identical |
| e_nonhuman | small tumor | a chrM target with `-only_target` (BT16) | Java and port: exit status 1 |
| e_threshold | small tumor | a threshold without `-write_off_target` (BT17) | Java and port: exit status 1 |
| e_nodir | small tumor | a missing output directory (BT19) | Java and port: exit status 0, errors logged, no files |
| small_t_th1 | small tumor | `-threads 1` | identical to the 2-thread run |
| wgs_t, wgs_n | a private whole-genome pair (80x tumor, 105x normal) | oncoanalyser command | run (section 3); Java references pending |

The synthetic BED is a regular grid on chr21 and chr22 plus a region nested in another (BT12), regions
1 and 3 bases apart (merged and kept apart) and pairs 50 and 60 bases apart (reads overlapping both,
BT18). Java's 2-thread partition_stats file has the rows out of partition order (chr1:1000000-1999999
at row 122 of the small tumor's), as BT9 says.

Corrected mode against Java (what the corrections change):

| Input | Files that change | What |
|---|---|---|
| small_t, small_n, e_partition | summary, coverage | BT1 (bases across partition starts), BT13 (UnmappedPercent 0.00334 to 0.00321 on small_t) |
| sbx_t, ultima_n | summary, coverage, flag_counts | as above; BT3: the unplaced unmapped reads move to QC-passed ("18576800 + 2381434 in total" becomes "20958234 + 0"; the mapped share 100.00% becomes 88.64%) |
| v_target | also frag_length, target_regions, off_target_overlaps, off_target_enriched | BT12 (the nested region keeps its 1500 bases: territory +900), BT5 (fragments whose lower read lies outside the targets: +10% fragments per grid region, 822 to 906; the off-target files get 7,999 overlap levels and 438 enriched sites (Java: headers only)) |
| v_target_only | also target_regions | BT12, BT5 through the upper read, BT18 (alone: 139 fewer reads, target rows unchanged) |
| e_noversion | gene_coverage, exon_coverage | BT10: GRCh38 gene coordinates, real depths in place of all zeros |
| e_nonhuman, e_threshold | | runs complete (exit 0) |
| wgs_t, wgs_n | summary, coverage, flag_counts, gene_coverage, exon_coverage | BT1, BT2, BT3, BT13 (section 2); frag_length unchanged |

Unit tests: 23 in bam-tools-rs (`cargo test --release`): the Java tests CoverageTest (filter types,
overlapping reads with mates and supplementaries), MetricsTest (read counts, dual-strand supplementaries)
and GeneCoverageTest (buckets, exon alignments) ported, plus tests of the findings, the record view, the
configuration, the off-target analyser (its fast overlap counts against Java's loop on 500 random
pile-ups) and the histogram. A small Java harness (compiled against the
jar) runs 704 synthetic cases through Java's CoverageMetrics (400 random depth distributions: mean,
median, MAD, standard deviation, filtered percentages and coverage levels, as raw double bits), the
coverage buckets (4 lengths) and OffTargetAnalyser (300 random pile-ups, which no BAM run reaches
because of BT5); `testdata/java_harness.tsv` holds them and two Rust tests replay them: all identical.

## 2. Findings

Corrected by default; `-java_compat` reproduces Java (except BT8 and BT9, which cannot be reproduced:
the port gives Java's single-thread result).

| Id | Java | Why it is wrong | Correction | Effect |
|---|---|---|---|---|
| BT1 | A read is counted only in the partition where it starts, and its bases past that partition's end are dropped (`processMatchedBases` breaks at the partition end; the next partition skips reads that start before it). | The first read length of every partition (3,102 partitions of 1 Mb for a genome) loses the depth of the reads crossing into it; the result depends on `-partition_size`. The same holds for mate-overlap removal across the boundary. | The reads that start before a partition and reach into it add their bases there (only per-base counts; the read-level counts stay with the partition the read starts in), with the fragment overlap removal applied to them too. | Small tumor: 536,013 more unfiltered bases (+0.008%). Whole genome: +16.8M (tumor) and +22.3M (normal) unfiltered bases, +0.007%; MeanCoverage 80.604 to 80.610 and 104.873 to 104.881, StdDevCoverage 36.084 to 36.081 and 25.031 to 25.021; DuplicatePercent and OverlappingReadPercent in the 5th decimal; the coverage distribution. |
| BT2 | Gene coverage matches a read against the exons overlapping the read's partition only. | A read starting near a partition's end that reaches into an exon starting in the next partition is not counted there (nor in the next partition, which skips it). | Every exon of the chromosome. | None on the small pair. Whole genome: 1 of 7,931 exons (one starting 81 bases past a partition boundary) gains 7% mean depth in both tumor and normal, and so its gene's row. |
| BT3 | The unplaced unmapped reads go to the flag counts' total with `passesQC = record.getReadFailsVendorQualityCheckFlag()`. | Inverted: reads that pass QC are counted as failed, and failed ones as passed. | `!fails`. | The total line and the mapped share. Whole genome: 54.5M (tumor) and 51.2M (normal) unmapped reads move from the QC-failed to the passed column; the mapped share goes from 99.57% to 96.88% (tumor) and from 99.57% to 97.54% (normal); BamFlagStats.mappedProportion, which ORANGE and qsee read, follows (0.9957 to 0.9688 and 0.9754). HG008 SBX tumor: 2.38M reads, 100.00% to 88.64%. The other lines are unchanged: Java adds the unmapped reads to the total only (whether they should count as primary, paired, read1/read2 as samtools does is unclear; kept). |
| BT4 | FlagStats: `isSecondary = read.getSupplementaryAlignmentFlag()`. | Copy error; for a consensus read the primary-duplicate decrement tests the supplementary flag twice. | The secondary flag. | None in practice (REDUX makes no secondary consensus reads). |
| BT5 | BamReader returns for every read outside the target regions before `checkTargetRegions`, the only caller of `OffTargetAnalyser.addRead`; reads that do reach it overlap a target by their own span, so `addRead` is never called. | `-write_off_target` writes empty files (header only): the off-target pile-up analysis never runs. Also a fragment whose first (lower) read lies outside the targets and whose mate overlaps one is never counted in the target region (the mate, being the upper read, is skipped as "already counted"). | Off-target reads go through `checkTargetRegions` as well (fragment counts for the targets they reach through the mate; the rest to the off-target analyser). With `-only_target`, where a lower read lying wholly before the region is never read for it, the upper read counts the fragment (its mate end from the MC tag). | Panels with `-regions_file`: off_target_overlaps and off_target_enriched get their data (synthetic targets on the small tumor: 7,999 overlap levels, 438 sites, against headers only); target_regions gains those fragments (+10%: 822 to 906 for a 1.5 kb region). The corrected `-only_target` and whole-partition modes then agree on the region counts where Java's two modes did not (247 and 294 for one region). |
| BT6 | OffTargetAnalyser processes a pile-up only when a later fragment starts past it, and nothing processes the last one of a partition. | Each partition's last pile-up is lost from the overlap counts and enriched sites. | The fragments still held are processed at the end of the partition. | With BT5 corrected. |
| BT7 | OffTargetAnalyser's peak region: when a new maximum depth is found the start moves but the end keeps the end of the previous (lower) plateau, or 0. | The peak can end before it starts (e.g. PeakPosStart 25, PeakPosEnd 24). | A new maximum sets both. | off_target_enriched rows. |
| BT8 | ExonCoverage depth arrays are shared by all partition threads and incremented without synchronisation. | A data race: an exon spanning a partition boundary can lose increments when the two partitions run at once; results can vary between runs. | Per-thread arrays, summed at the end. | Not observed in the runs so far (needs simultaneous increments of one base); not reproducible, so Java mode gives the race-free result. |
| BT9 | partition_stats, target_regions and off_target_enriched rows are written as threads finish partitions. | Row order varies between runs with more than one thread. | Partition order (Java's order with one thread). | Row order only. |
| BT10 | Without `-ref_genome_version`, MetricsConfig takes the version from the BAM header, but the Ensembl cache takes `RefGenomeVersion.from(configBuilder)`, which defaults to V37. | A GRCh38 BAM gets GRCh37 gene coordinates (or the unversioned files read as 37, with chromosomes without "chr"), so gene and exon coverage come out empty (all zero depth). | The BAM's version. | Only without `-ref_genome_version` (oncoanalyser always passes it). |
| BT11 | BaseCoverage and the read-group coordinates handle M, I, D, N, S, H; `=` and `X` fall to the default case and consume nothing. | For an aligner writing `=`/`X` the later blocks are placed at wrong positions and those bases are not counted. | `=` and `X` count as M. | None for BWA (M only). |
| BT12 | hmf-common BaseRegion.checkMergeOverlaps (used on `-regions_file` regions) sets the merged end to the next region's end even when it is smaller. | A region lying inside the one before cuts the outer region short (100-500 and 200-300 become 100-300). | The larger end is kept (as ChrBaseRegion.checkMergeOverlaps does). | Target-region files with nested regions. |
| BT13 | The MATE_UNMAPPED filter counts the read's length (soft clips included); the other read-level filters (low map quality, duplicates) count aligned (M) bases. | Mixed units in the summary's filtered percentages. | Aligned bases. | MATE_UNMAPPED bases -4%: UnmappedPercent 0.00334 to 0.00321 (small tumor), 0.00337 to 0.00324 and 0.00335 to 0.00329 (whole-genome tumor and normal). |
| BT14 | The median (and the MAD) take the value at rank total / 2 of the depth distribution. | For an odd count this is one below the median, and for one base it is 0 whatever its depth. | Rank (total + 1) / 2. | Only visible on tiny regions. |
| BT15 | The collapsed coverage distribution (max coverage above 3000, the panel default of 10000) puts coverage 1 in the first bucket, labelled 2, because `getNextCoverageBucket` returns the first bucket for levels below it but the bucket below for every other level. | No row for 1x, and the 1x bases are counted under 2x. | Level 1 has its own row. | Panels (coverage file with `-regions_file`). |
| BT16 | With `-only_target`, a target region on a non-human contig (chrM, an alt contig) makes ChrBaseRegion.overlaps call HumanChromosome.fromString, which throws; the partition thread exits the program (status 1). | A crash on valid input. | Non-human contigs compare by name. | Panels with such regions. |
| BT17 | With `-regions_file` and `-off_target_frag_overlap_threshold` > 0 but no `-write_off_target`, the enriched-region file is opened but there is no OffTargetAnalyser: a NullPointerException after the first partition (status 1). | A crash on a valid combination of options. | No analyser, no rows; the file has its header. | That option combination. |
| BT18 | With `-only_target`, each target region reads every read overlapping it; a read overlapping two regions (closer than a read length) is counted in both: total and duplicate reads, flag counts, fragment lengths, read-level filter counts and gene coverage. | Double counting. | The read-level counts in the first region it overlaps; its bases and the per-region fragment counts still count in each region. | Synthetic targets on the small tumor (two pairs 50 and 60 bases apart): 139 fewer reads; target_regions unchanged. |
| BT19 | BamMetrics does not create the output directory; without it every writer logs an error and the run ends with status 0 and no files. | Silent failure. | The directory is created. | Only for a missing directory (oncoanalyser creates it). |
| BT20 | hmf-common DriverGeneFile.write puts reportLowExpression's value under the reportHighExpression header and the other way round (header in enum order, values in toString order). | A panel written by Java reads back with the two swapped. | hmf-common-rs `write_driver_genes` writes each value under its own column; `write_driver_genes_java` reproduces Java. | Not used by BamTools (shared code for later ports). |

Remarks, not findings:

- `findOverlappingFragmentStart` tests `suppData.Position >= suppData.Position` (always true) where
  `>= readStart` was meant; `positionWithin` checks the same, so there is no effect.
- `checkMergeOverlaps` merges regions with a one-base gap (end >= next start - 2); kept.
- `FragmentLengths.mMaxReadLength` is computed and never used.
- The unmapped-read query opens the BAM with htsjdk's default STRICT validation (the partition readers
  use SILENT), so an invalid unmapped record would end the Java run; htslib does not validate.
- OffTargetAnalyser is quadratic in a pile-up (a name scan per read, pairwise overlaps, a depth sum over
  all fragments for each base of the peak). On a panel's sparse off-target reads that is cheap, but on
  continuous coverage a whole partition is one pile-up: the first corrected run with targets on WGS-like
  data had not finished after 10 minutes. The port computes the same values in O(n log n) (a name set,
  binary searches on the sorted starts and ends, a difference array), checked against Java's loops on
  the 300 harness cases and 500 random pile-ups.
- With `-only_target` a region reads only the reads overlapping it, so a fragment spanning a region with
  neither read inside it is not counted there (it is without `-only_target`): 170 against 179 fragments
  for a 100-base synthetic region. A design limit of the mode, not corrected.
- Java's two modes give different target region counts for the same region (BT5, BT18); corrected,
  they agree except for the fragments of the remark above.

## 3. Timings and memory

One machine (Apple M-series, 20 cores, 128 GB, shared: other ports' runs were going on), inputs on a
USB SSD, `/usr/bin/time -l`, 2 threads for both, Java `-Xmx16g`. CPU = user + sys. All in Java mode.

| Run | Java wall | Java CPU | Java peak | Rust wall | Rust CPU | Rust peak | Wall speed-up |
|---|---|---|---|---|---|---|---|
| small_t | 105.8 s | 190.6 s | 2.23 GB | 12.5 s | 23.3 s | 77 MB | 8.4x |
| small_n | 118.9 s | 213.4 s | 2.25 GB | 14.5 s | 26.8 s | 77 MB | 8.1x |
| sbx_t | 55.3 s | 77.5 s | 2.68 GB | 9.1 s | 13.7 s | 79 MB | 6.1x |
| sbx_n | 25.1 s | 47.1 s | 2.84 GB | 3.8 s | 6.5 s | 79 MB | 6.6x |
| sbx_n_cram | 43.6 s | 85.9 s | 2.98 GB | 9.5 s | 16.9 s | 118 MB | 4.6x |
| ultima_n | 119.8 s | 165.4 s | 3.18 GB | 19.7 s | 32.7 s | 79 MB | 6.0x |
| ultima_t | 59.9 s | 101.5 s | 3.05 GB | 11.5 s | 18.4 s | 79 MB | 5.1x |
| v_quals | 47.9 s | 90.4 s | 2.32 GB | 6.3 s | 11.8 s | 71 MB | 7.6x |
| v_nodup_zero | 36.8 s | 71.7 s | 2.33 GB | 5.6 s | 10.1 s | 71 MB | 6.5x |
| v_regions | 13.8 s | 30.1 s | 2.28 GB | 2.2 s | 3.8 s | 71 MB | 6.2x |
| v_target_only | 2.6 s | 7.8 s | 1.39 GB | 1.0 s | 1.0 s | 71 MB | 2.6x |
| v_target | 32.6 s | 68.8 s | 2.35 GB | 8.4 s | 16.2 s | 76 MB | 3.8x |
| v_threads1 (1 thread) | 76.0 s | 80.9 s | 2.32 GB | 12.0 s | 12.0 s | 79 MB | 6.3x |
| e_noversion | 42.8 s | 88.7 s | 2.30 GB | 6.3 s | 11.8 s | 73 MB | 6.7x |
| e_partition | 34.8 s | 72.8 s | 2.33 GB | 5.6 s | 10.5 s | 71 MB | 6.1x |
| wgs_t (80x tumor, 2.0G reads) | pending | | | 432.2 s | 777.1 s | 109 MB | |
| wgs_n (105x normal, 2.5G reads) | pending | | | 548.9 s | 970.4 s | 105 MB | |

The whole-genome runs are the final binary's; earlier runs (unmapped reads read last, on one thread) took 453.5 s and 554.6 s, and two runs of each BAM gave identical files. The corrected mode costs about the same (wgs_t 467.5 s, small_n 14.8 s). Threads: the small tumor takes
24.3 s at 1 thread and 12.5 s at 2 (1.9x). The pipeline figure for Java (about 8 CPU hours for a 100x
pair, `../RUST_PORTS.md`) against the port's 0.49 CPU hours for this pair is a factor of about 17 in CPU;
the Java references will give the measured one.

## 4. Speed: what was done and where the time goes

- Workers take partitions in order from a shared counter, each with its own indexed reader and BGZF
  decompression (libdeflate); no per-read allocation (the record is reused; the read view borrows its
  data; one pass over the aux data, memchr for strings, finds CR, UT, SA and MC; read-group coordinate
  vectors are recycled); per-thread accumulators, summed at the end.
- htslib's default 32 KB file reads left both threads waiting on the USB volume: 2 threads ran no faster
  than 1 (31.7 s against 32.6 s on the small tumor). 1 MB reads (`-file_buffer_kb`, default 1024, as
  sage-rs) fixed it (4 MB gave 15.9 s, 1 MB 18.5 s on a busy machine; the later small_t run 12.5 s).
- The per-base loop is branch-free (quality test, cap test, increment) and the coverage histogram of
  every partition base runs on four sub-histograms with whole zero blocks skipped; before that change
  the histogram alone (3.1 billion bases) cost more than decompression (1 thread: 32.6 s to 24.3 s).
- Exons ending before a read are skipped with a running maximum of exon ends (Java scans the partition's
  exon list from its start for every read).
- The Ensembl transcript file is streamed: peak memory went from 500 MB to 77 MB.
- The off-target analyser is O(n log n) (section 2, remarks).
- The unplaced unmapped reads (65M in the whole-genome tumor: 22 s on one thread after all partitions,
  as Java does it) are a task of their own, taken first by one worker while the other starts on chr1
  (their counts are sums, so the order changes nothing).
- Profile of the whole-genome tumor (macOS `sample`, 15 s of both threads): libdeflate decompression
  with its Huffman table builds and CRC about 32%, the `read` system call 9% (the BAM streams from
  the USB volume at about 210 MB/s), the per-base counting 7%, the record view 5%, htslib record reading
  and copying 9%, read groups 4%, flag statistics 2%; the rest is spread thin. The small tumor at 1
  thread: decompression 45%, per-base counting 10%, record view 10%, htslib reading 15%.
- What is left is mostly decompression, which any whole-BAM reader pays; sharing one pass over the BAM
  with AMBER, COBALT and TEAL (`../RUST_PORTS.md`) is where the next factor would come from.

## 5. Open points

- The whole-genome comparison with Java, once the Java references of the two BAMs exist (the port's
  outputs in both modes are made). The partition_stats rows of a 2-thread Java run are compared sorted
  and without ProcessTime.
- `EnrichedRegionRepeatMasks` (a separate main class, with hmf-common RepeatMaskAnnotations) is not
  ported; neither BamMetrics nor oncoanalyser calls it.
