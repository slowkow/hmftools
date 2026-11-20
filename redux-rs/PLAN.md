# redux-rs: status

A port of REDUX 2.1 (`../redux`). Parity means: the BAM holds the same records in the same order as
Java's, and the other output files are byte-identical. See README.md for how the thread count enters.

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
| Unit tests | the Java tests of the ported classes (71 core, 36 SBX, 21 UMI, BQR, jitter, Ultima); 184 Rust tests | `tests/`, and next to the modules |

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
- A whole-genome run from raw alignments, a Linux build, reference genome 37.

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

## Measurements

Quiet machine (20 cores), 16 threads:

| Input | Java wall / CPU / memory | Rust wall / CPU / memory |
|---|---|---|
| Illumina raw, chr21 + chr22 + unplaced reads, 111 M reads (40 M unmapped in the first pass) | 879 s / 1,389 s / 23.5 GB | 256 s / 603 s / 0.5 GB |
| Illumina, two chromosomes of a whole genome, 343 M reads, 15 GB (Java run under a profiler) | 1,123 s / 4,404 s / 23.5 GB | 525 s / 2,073 s / 0.5 GB |

Both inputs fill only a few of the 32 partitions, so neither tool uses 16 cores here. A whole genome
from raw alignments has not been run.

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
