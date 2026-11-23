# redux-rs: status

A port of REDUX 2.1 (`../redux`). By default the port corrects the faults found in the Java code
(REPORT.md, section 3); `-java_compat` reproduces them. Parity means, with `-java_compat`: the BAM holds
the same records in the same order as Java's, and the other output files are byte-identical. See
README.md for how the thread count enters. Every comparison with Java in this file was made in Java
mode (before the corrections existed, the only mode).

## Done and compared with Java

| Part | Java | Rust |
|---|---|---|
| Command line, partition split, read-length sampling | ReduxConfig, ReduxApplication, PartitionThread | `config.rs`, `main.rs`, `regions.rs` |
| Input, several BAMs merged by position | BamReader | `bamin.rs` |
| Unmapping in problematic regions, both passes | ReadUnmapper, RegionUnmapper | `unmap.rs`, `stage1.rs` |
| Duplicate groups | FragmentCoords, ReadCache, DuplicateGroup, Statistics | `dup.rs` |
| Consensus reads | ConsensusReads, BaseBuilder, IndelConsensusReads, ReadParseState, IlluminaRoutines, TemplateReads | `consensus.rs` |
| Sorted output, concatenation, index | SortedBamWriter, FileWriterCache, samtools | `sorted.rs`, `out.rs` |
| Base quality recalibration | package `bqr` | `bqr.rs` |
| Jitter analysis, MSI prediction | packages `jitter`, `ms_model` (not the training tool) | `jitter/` |
| SBX | SbxRoutines, SbxDuplexIndelBuilder, SbxIndelConsensus, SbxDuplicateCollapser, SbxBamSampler | `sbx.rs`, `sbx_dup.rs`, `consensus/sbx_indel.rs` |
| Ultima | UltimaRoutines, UltimaStats | `ultima.rs` |
| UMIs | UmiConfig, UmiGroupBuilder, PolyGUmiMerger, JitterUmiMerger, UmiStatistics | `umi.rs` |
| Unit tests | the Java tests of the ported classes (71 core, 36 SBX, 21 UMI, BQR, jitter, Ultima) and tests of the corrected paths; 211 Rust tests | `tests/`, and next to the modules |
| The switch between corrected and Java behaviour | not in Java | `compat.rs` |

## Corrections

Done: 24 of the 40 numbered findings are corrected by default, each behind `compat::java("<id>")` with
the Java path beside it (REPORT.md, sections 3.1 to 3.5). The Java-mode regression (88 stored Java
baselines) gives the same result as before the corrections.

Left as in Java, in both modes, with the reason in REPORT.md section 3.6:
- R12 (the outlier pass of the jitter tables) and two of the three options of R13: the code holds no
  rule to apply.
- R15 (duplicates are not matched across a partition boundary): needs the partitions to read and
  group across their boundaries; the redesign is described in section 3.6.
- R23 (Ultima consensus on bases that disagree): Java has no rule; the port stops with a message.
- R30 (`-output_id` and the names of the side files): downstream tools find these files by name.

## Not done

- UMIs on real data: compared with Java on simulated UMIs only. UMIs with SBX or Ultima reads are refused.
- CRAM input: htslib reads it, but the result has not been compared with Java, whose CRAM decoder treats
  some attributes differently.
- `-read_output` (the per-read TSV), `-run_checks` logging, `-log_read_ids`, `-perf_log_time`,
  `-log_dup_group_size`, `-log_read_cache`: accepted, no effect, apart from the SBX consensus fallback
  under `-run_checks`.
- The R plots (`-bqr_write_plot`, `-write_msi_plots`): accepted with a warning.
- The separate tools in the same jar: the microsatellite site finder (`ms_sites`) and the MSI model
  training (`MsModelBuilder`).
- A whole-genome run from raw alignments, reference genome 37.

## Design notes

- **Records stay htslib records.** A read is decoded once by htslib, edited in place and written back.
  Fields are read from the raw record; nothing is copied into a second read type.
- **htsjdk's write rule is reproduced.** htsjdk writes a record's stored bytes unchanged until its name,
  bases, qualities, cigar or an attribute is set, and from then on writes the attributes sorted by tag
  with integers in their smallest type. `rec.rs` does the same (`set_attribute`, `decode_attributes`).
- **java.util.HashMap order is reproduced** where Java iterates one and the order reaches the output:
  the read cache's fragment map, the Ultima statistics.
- **First pass.** Reads unmapped in the problematic regions go to one file per partition, tagged with
  their order of discovery; each file is sorted (four at a time, in memory) and then streamed into its
  partition as one more input. Java writes one BAM, sorts it with samtools and reads it back.
- **Output.** Each partition writes a headerless block-compressed record stream; the final BAM is the
  header plus those streams copied block by block. Records are reused after writing.
- **Jitter.** Sites are sorted per chromosome and each worker follows the sorted reads with a cursor;
  counts go to a dense per-worker window. Java queries an interval tree and allocates per read.
- **BQR.** A flat slot array per position instead of a map of lists per position; per-reader key table.
  In Java mode the array is Java's block of 100K positions. In the default mode it is a window that
  moves with the reads (R2), and the positions next to a partition boundary are set aside by both
  partitions and counted by whichever completes second (R16), so no partition waits for another.
- **Read cache.** In Java mode a partition whose cache nears 100,000 reads waits for the read counts of
  the earlier partitions, to give Java's one-thread result. In the default mode the counter starts
  again with each region (R21) and nothing waits.
- **The switch.** `compat::java("<id>")` is an atomic load in the default mode and in Java mode (a set
  lookup only under `-java_compat_except`). The BQR reader, the read cache and the partition reader
  read it once when they are made.

## Measurements

Quiet machine (20 cores), 16 threads:

| Input | Java wall / CPU / memory | Rust wall / CPU / memory |
|---|---|---|
| Illumina raw, chr21 + chr22 + unplaced reads, 111 M reads (40 M unmapped in the first pass) | 879 s / 1,389 s / 23.5 GB | 256 s / 603 s / 0.5 GB |
| Illumina, two chromosomes of a whole genome, 343 M reads, 15 GB (Java run under a profiler) | 1,123 s / 4,404 s / 23.5 GB | 525 s / 2,073 s / 0.5 GB |

Both inputs fill only a few of the 32 partitions, so neither tool uses 16 cores here.

Whole genome: a vendor (DRAGEN) normal BAM, 86 GB, 1.56 billion reads, read from and written to a
network disk, so the 16-thread wall times show the disk, not the tools:

| Threads | Java wall / CPU / memory | Rust wall / CPU / memory |
|---|---|---|
| 16 | 2,215 s / 13,488 s / 35.0 GB | 2,058 s / 10,909 s / 0.9 GB |
| 1 (both tools running at the same time) | 15,122 s / 12,919 s / 39.6 GB | 9,906 s / 7,792 s / 0.8 GB |

At 16 threads both tools spent the same time in each phase (first pass about 9 minutes, partitions 16
to 18, joining and indexing 9), all bound by the disk. The 16-thread Rust run was made before the index
was built while writing. The output BAM was 77.7 GB against Java's 77.4 GB. The Linux runs below
are on a local disk.

After the second optimisation pass, the same whole genome on macOS at 16 threads from an external SSD
took 1,013 s (5,469 s CPU, 1.1 GB): first pass 2.7 minutes, partitions 10, joining the partition files
and the rest 4.2. All 1,780,126,609 records and the four side files are identical to Java's one-thread
output. The Java and Linux figures in this section were measured before that pass.

Since then the partition files are appended to the output while later partitions run (and deleted as
they are appended). On the same run this saved no time (1,034 s; identical output): the external
disk writes about 300 MB/s of compressed data, which bounds the second copy of the output whenever
it is made. A third pass over the profile took 6% of the instructions off the one-thread test.

Linux (x86-64, 16 virtual CPUs, 128 GB, local SSD; built with the stock Debian 12 toolchain), same
vendor data, Java 17:

| Input | Threads | Java wall / CPU / memory | Rust wall / CPU / memory |
|---|---|---|---|
| chr21 slice | 4 | 385 s / 429 s / 33.8 GB | 153 s / 154 s / 0.4 GB |
| chr21 + chr22 + unplaced reads | 16 | 1,295 s / 2,131 s / 33.2 GB | 401 s / 929 s / 0.9 GB |
| Whole genome, 86 GB | 16 | 3,030 s / 31,074 s / 35.7 GB | 1,583 s / 16,693 s / 1.6 GB |

CPU is user time. In the whole-genome run Rust's phases were 3 minutes for the first pass, 20 for the
partitions and 4 for joining, jitter and the index; Java's were 5, 32 and 8. The BAMs match except
for the order of unplaced reads and, in the whole genome, the same few hundred duplicate flags on one
very deep contig that differ on macOS at 16 threads (Java's result there depends on thread timing).

Both modes on the same whole genome (2026-11-22, shared machine, 32 partitions, 6 worker threads,
external SSD): REPORT.md, section 6.

One thread, taken while other jobs were running; CPU time is the fairer figure:

| Input | Java wall / user / memory | Rust wall / user / memory |
|---|---|---|
| Illumina, raw bwa output, 25.8 M reads | 308 s / 319 s / 20.1 GB | 139 s / 135 s / 0.5 GB |
| The same with duplex UMIs (simulated) | 1,269 s / 510 s | 199 s / 175 s / 0.5 GB |
| SBX tumor chr21, 23.8 M records out | 470 s / 488 s / 18.6 GB | 240 s / 225 s / 0.6 GB |
| Ultima tumor chr21, no consensus | 1,453 s / 727 s | 261 s / 233 s / 0.4 GB |
| BQR and jitter only, Illumina slice | 194 s | 24 s |

Where the Rust time goes on Illumina (one thread): output compression about a third, BQR 11%,
allocation 10% before records were reused. Output compression level, same input:

| Level | Wall | BAM size |
|---|---|---|
| 5 | 148 s | 1,148 MB |
| 4 (default) | 138 s | 1,178 MB |
| 3 | 123 s | 1,236 MB |
| 2 | 120 s | 1,267 MB |
| Java (zlib 5) | | 1,181 MB |
