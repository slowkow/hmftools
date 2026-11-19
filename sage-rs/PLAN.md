# SAGE reimplementation plan

Target: a native reimplementation of SAGE 5.1 (`sage/`, upstream `a224a2c`) that produces the same
calls as the Java version on tumor/normal Illumina data, and runs substantially faster.

## Status

The Rust port runs the full tumor/normal Illumina somatic pipeline and writes a bgzipped, tabix-indexed
VCF. On every region tested so far the decompressed VCF is byte-identical to Java SAGE 5.1 run with
`-threads 1`, at any Rust thread count.

| Step | State |
|---|---|
| 0. Baseline and harness | Done: `bench/run_java.sh`, `bench/run_rust.sh` (same arguments, oncoanalyser somatic settings minus TINC). Comparison is a plain `diff` of the decompressed VCFs. |
| 1. Skeleton | Done: CLI, reference data, partitioning, BAM slicing, VCF writer. |
| 2. Candidate stage | Done. |
| 3. Evidence stage | Done, including fragment sync, BQR and the MSI jitter model. |
| 4. Filters, phasing, dedup | Done. |
| 5. Parity sweep | Whole genome compared with and without TINC (below), plus chr21, HLA region, chrM and part of chrX against single-threaded Java. chr21 is also byte-identical in `-high_depth_mode`, `-germline`, tumor-only and `-panel_only` runs. The Java unit tests are ported (`src/tests`), with further tests for the port's own alignment reading and fast paths; 157 tests pass. |
| 6. Optimisation | Multi-threading; table-based base decoding; the BAM tags are parsed in place; read buffers are reused in both stages; the tumor's reads are read once per partition and kept for the evidence stage; overlapping plain-match mates are merged by a fast path; the candidate stage skips from mismatch to mismatch. About half the CPU of the first complete version (measurements below). |
| 7. TINC | Done: `-run_tinc` and the standalone `sage-rs tinc`. |
| 8. SBX | Done: `-sequencing_type SBX`. chr21 of the public HG008 tumor/normal pair (119x / 44x, after Java REDUX) is byte-identical to Java, with and without soft filters (3,316 and 68,809 records). |
| 9. Ultima | Done: `-sequencing_type ULTIMA`. chr21 of the public HG008 ppmSeq tumor/normal pair (157x / 160x, after Java REDUX without `-form_consensus`) is byte-identical to Java, with and without soft filters (2,689 and 87,238 records). |
| 10. Append, split reads, MSI samples, multi-sample | Done: `sage-rs append`, `SplitReadSegment`, MSI-specific jitter parameters, several tumor / reference samples. Compared with Java on the public HG008 data (chr21): append onto SBX and Ultima VCFs with and without soft filters, two samples at once, `-high_depth_mode` and an older-format VCF; simulated MSI predictions for SBX, Ultima and the default path; two tumors and/or two references. All byte-identical. Split reads have unit tests only (no RNA BAM compared). |
| 11. CRAM | Done: CRAM 3.0 and 3.1 input, decoded against `-ref_genome`. chr21 from CRAM is byte-identical to Java reading the same CRAM for the Illumina, SBX and Ultima pairs (5,908, 3,316 and 2,689 records). For the Illumina pair that is not the BAM result, in Java either: see QUIRKS.md. |
| 12. Smaller options | Done: `-write_frag_lengths` and `-specific_var_file`, in the caller and in append, and `-run_tinc` with two tumors. Compared with Java on the public HG008 data and the Illumina pair: the fragment length files (1.1 to 5.3 million lines) and the VCFs are byte-identical. One site in a two-tumor `-disable_soft_filter` run differs because Java's own answer there is arbitrary (QUIRKS.md, `IndelDeduper`). |
| 13. Visualiser | Done: `-vis_variants`, `-vis_pass_only`, `-vis_max_support_reads`, `-vis_output_dir`, `-vis_purple_vcf` (with the amino acid track). Compared with Java on the public HG008 data (SBX and Ultima, one and several samples, tumor only) and on the Illumina pair (paired reads, merged fragments, a PURPLE VCF): the same file names, and the same bytes once the read rows that Java orders at random are sorted. The VCF of a run with the visualiser on is identical to Java's. 181 unit tests pass. See "Visualiser" in QUIRKS.md for what cannot match. |

Regions compared (records identical / total), tumor and normal WGS BAMs on an external volume, 20-core machine:

| Input (one thread) | Records | Java | Rust |
|---|---|---|---|
| Illumina chr21 | 6,023 / 6,023 | 501 s, 3.6 GB | 59 s, 150 MB |
| Illumina chr21, CRAM | 5,908 / 5,908 | 767 s | 110 s |
| SBX chr21 | 3,316 / 3,316 | 106 s | 37 s |
| SBX chr21, CRAM | 3,316 / 3,316 | 310 s | 81 s |
| Ultima chr21 | 2,689 / 2,689 | 253 s | 109 s |
| Ultima chr21, CRAM | 2,689 / 2,689 | 692 s | 201 s |
| Illumina chrM | 128 / 128 | 25 s | 4 s |

Before the second optimisation round the Rust times were 114 s (Illumina), 61 s (SBX) and 172 s
(Ultima). Other Illumina regions compared with the earlier build: chr21 with `-disable_soft_filter`
(110,762 records), the HLA region of chr6 (851) and chrX:1-6 Mb (482).

chr21 with more threads (first complete version): Java 8 threads 72 s (8.5 GB); Rust 8 threads 19 s (0.5 GB); Rust 16 threads 12 s.
Java's own 8-thread output differs from its 1-thread output on 152 lines (phase set ids, as section 4
predicts); the Rust output is the same at 1, 8 and 16 threads.

Whole genome, 16 threads each (Java with `-Xmx24g`): Java 2,450 s and 26.9 GB peak; Rust 812 s and 0.9 GB.
Both write the same 360,440 records. 355,867 lines are identical; the other 4,573 differ only in the
`LPS` id numbers, which is the Java thread-order effect (`bench/compare_vcf.py`). No other field differs.

With `-run_tinc` (the oncoanalyser setting), whole genome, 16 threads each: Java 3,407 s, Rust 764 s
(both measured while other jobs were running, so only the ratio is meaningful). The pre-TINC VCFs have
the same 611,876 records and the TINC-rewritten VCFs the same 365,286; in both, the only differing
lines are the 4,704 with thread-order `LPS` ids. The sample's TINC level is 0, so recovery was checked
separately: `bench/simulate_tinc.py` adds simulated contamination to a VCF, and the standalone Java and
Rust analysers give byte-identical output at 2%, 10% and 50% (5,132 to 14,381 variants recovered).
Standalone TINC on the whole-genome VCF takes 113 s and 17 GB in Java, 14 s and 0.5 GB in Rust.

After the second optimisation round, the same `-run_tinc` whole genome at 16 threads takes 593 s with a
1.7 GB peak and writes the same two VCFs as the build compared with Java above. CPU time fell from
9,030 to 5,380 core-seconds; the wall time fell less (771 s before) because the run now waits on the
disk: it averages 9 of 16 cores, and the same build runs four chromosomes at 16 cores once the BAMs are
in the page cache. `-perf_summary` reports how long the workers waited on the writer (54 s summed over
16 threads here) and how long the writer was busy (64 s). The peak is above the earlier 0.9 GB because
of the read cache (32 MB per thread by default; `-read_cache_mb 0` turns it off).

Rust CPU profile on Illumina chr21 after the second round (one thread): BAM decompression and reading
34%, matching reads to candidates and quality scores 13%, allocation 11%, decoding records 8%, the
jitter check 7%, read-context building 6%, fragment sync 5%. For Ultima, read-context building (core
extension over homopolymers, repeat search) is about half.

Tried and dropped: the mimalloc allocator was 10 to 15% faster, but its memory grew throughout a
whole-genome run (12 GB against 2 GB), because region results are allocated by worker threads and
freed by the writer thread.

Java CPU by stage on chr21 (flight recorder): candidate stage 63%, evidence stage 34%, everything after
that under 1%. BAM decoding inside htsjdk is about 20% of the total.

## 1. What SAGE does

Per chromosome, the genome is cut into 100 kb partitions (`pipeline/ChromosomePartition`). A pool of
threads pulls partitions from a queue; each `RegionTask` runs:

1. **Candidates** (`candidate/`) - slice the tumor BAM, walk each read's CIGAR against the reference,
   record every SNV / MNV (up to 3 bases) / indel / soft-clip-implied insert in a position-indexed
   evicting array (`RefContext` -> `AltContext`), build a read context (core + flanks, repeats,
   microhomology) for each alt, keep those with enough support or at hotspots.
2. **Evidence** (`evidence/`, `common/ReadContextMatcher`, `quality/`, `sync/`) - slice the tumor BAM
   a second time, then the reference BAM. Overlapping mates are merged into one synthetic read
   (`sync/`). Each read is tested against every candidate it covers: locate the variant index in the
   read, compute base/map quality (BQR table, MSI jitter model, edge and event penalties), classify
   as full / partial-core / core / realigned / ref / jitter, and accumulate counters, strand bias,
   edge distance, fragment coordinates, and phasing groups.
3. **Filter, phase, dedup** (`filter/`, `phase/`, `dedup/`) - hard filters, tiered soft filters
   (binomial tests), local phase sets, MNV/SNV/indel de-duplication.
4. **Write** (`vcf/`) - a reorder buffer emits regions in order to a bgzipped, indexed VCF.

Size: 26.7k lines of Java in `sage/src/main`, 9.1k lines of tests.

| Package | Lines | First milestone? |
|---|---|---|
| candidate, evidence, common, pipeline, select | ~8.7k | yes - the hot path |
| quality, sync | ~2.1k | yes |
| filter, phase, dedup, vcf | ~5.2k | yes |
| top-level config and reference data | ~1.3k | yes |
| seqtech (Illumina part) | ~0.5k | yes |
| seqtech (Ultima, SBX) | ~2.9k | no - reject `sequencing_type` other than ILLUMINA |
| tinc | ~1.4k | second milestone (needed if the pipeline passes `run_tinc`) |
| append | ~0.6k | later - reuses the evidence stage |
| vis, utils | ~3.8k | no (vis was added last: status row 13) |

So the first milestone is about 18k lines of Java plus a small slice of `hmf-common`
(`Doubles`, `BaseQualAdjustment`, BQR / jitter file readers, BED / hotspot readers, HLA gene
coordinates, `SamRecordUtils`) and replacements for htsjdk, commons-math3 and Guava.

## 2. Language: Rust

- htslib bindings give BAM, CRAM, indexed FASTA, bgzf and tabix with libdeflate, which is the same
  I/O layer a C++ port would use.
- The concurrency model (independent region tasks, one ordered writer, one shared phase-set counter)
  maps directly onto a work-stealing pool, and the compiler rules out the data races the Java code
  currently has (section 4).
- `cargo` 1.97 is already installed here; the Java unit tests port naturally to `#[test]`.

C++ would reach the same speed. It offers no advantage for this codebase and costs more in build
setup and memory-safety review.

## 3. Where the speed should come from

From reading the code (the profile in the Status section agrees on where the time goes):

- **The tumor BAM is decompressed and parsed twice** (candidate stage, then evidence stage), through
  htsjdk, which inflates BGZF single-threaded per reader and materialises a `SAMRecord` per read.
- **Per-read allocation.** Candidate stage creates `String`s and an `AltRead` for every mismatch
  (`RefContextConsumer.processAlignment`). Fragment sync builds a whole new `SAMRecord` through
  `SAMRecordSetBuilder` and copies every tag for each overlapping pair (`sync/CombinedSyncData`).
- **Per read x variant overhead.** BQR lookup concatenates a `String` key and probes a `HashMap`
  (`ReadContextQualCache.getQual`); the consensus-type tag is re-read on each call; HLA overlap is
  re-evaluated per read although it depends only on variant position.
- **Phasing** does linear scans with `List.contains` per read and quadratic merges per region.
- **Scheduling.** Threads are created per chromosome and each opens its own BAM readers; partition
  size is tied to `min(threads, 4)`.

A straight port removes most of the allocation and lookup cost. Reading the tumor BAM once per
partition and sharing decoded records between the two stages is the main structural change, and
comes only after parity is established.

Measured results are in the Status section. The port so far is a straight translation with
threads; none of the structural changes above have been made yet.

## 4. What "same output" means

Byte-identical VCFs are not a usable target, because the Java version is not reproducible against
itself when `threads > 1`:

- `LPS` ids come from one counter shared by all threads (`phase/PhaseSetCounter`), so they depend on
  scheduling.
- `NEAR_HOTSPOT` is set after the variants are handed to the writer (`RegionTask.finaliseResults`),
  so it races with output.
- `IndelDeduper` orders tied combinations by identity hash, which varies between JVM runs.

Parity definition, against Java run with `threads 1`:

1. Same set of records (CHROM, POS, REF, ALT) and same FILTER per record.
2. Integer INFO / FORMAT fields equal.
3. Float fields equal as printed; any difference must be traced to a last-digit libm or formatting
   difference and listed.
4. `LPS` equal up to renumbering.

Things that need care to meet that:

- Java `Math.round` is `floor(x + 0.5)`; Rust `round` differs on negative halves (affects QUAL).
- Java `Math.min` / `max` propagate NaN; several filter formulas divide by zero on purpose.
- `(int)` casts truncate and saturate; `int *= double` truncates (`ReadContextCounter.tumorQuality`).
- Binomial CDF / PMF must follow commons-math3 3.6.1 (`Beta.regularizedBeta`, saddle-point PMF).
- Sorts must be stable. Three places depend on Java `HashMap` iteration order
  (`PhasingUtils.removeUninformativeLps`, `FragmentSync.emptyCachedReads`, `AppendVariantPhaser`).
- htsjdk VCF number formatting and INFO / FORMAT key ordering.
- Read order and filtering from `BamSlicer` (`queryOverlapping`, drops duplicates and secondaries,
  keeps supplementaries only in the candidate stage), and the max-depth caps that depend on it.

Known oddities in the Java code are ported as they are, and listed in `QUIRKS.md` as they are found,
so they can be raised upstream rather than silently fixed. Examples so far: lost phase groups on key
collision in `PhasedGroupCollection.merge`; an always-true check in `PhasingUtils.mergeUninformative`;
`RC_JIT` declared with 3 values but written with 2; synced reads dropped on exception.

## 5. Steps

**0. Baseline and harness** (before any Rust)
- Build the Java jar (`mvn -pl sage -am package -DskipTests`).
- Pick test data: one tumor/normal pair with REDUX outputs (`*.redux.bqr.tsv`,
  `*.redux.jitter_params.tsv`, `*.redux.ms_table.tsv.gz`), plus the hotspot, panel and
  high-confidence resources. Take the exact command line from the oncoanalyser SAGE module.
- Record golden VCFs with `threads 1`: a handful of small regions (dense panel gene, HLA,
  a microsatellite-rich region, chrM) and one whole chromosome.
- Profile the Java run (JFR) and record wall time and CPU by stage. This decides the order of
  optimisation work in step 6.
- Write the VCF comparison tool implementing section 4.

**1. Skeleton** - CLI with the same option names and defaults; reference data loading (FASTA, BED,
hotspots, BQR, jitter files, panel merge); partitioning; BAM slicing with the same filters; ordered
VCF writer. Check: header matches, empty body.

**2. Candidate stage** - CIGAR walk, evicting array, read-context builder, repeats, microhomology,
tiering, depth caps. Check against Java `log_candidates` output, which isolates this stage.

**3. Evidence stage** - matcher, raw context, quality calculator, jitter, realignment, fragment
sync, counters. This is the largest and most intricate step. Check against the full VCF with
`disable_soft_filter`.

**4. Filters, phasing, dedup** - including the commons-math3 numerics. Check against the full VCF.

**5. Parity sweep** - whole genome, tumor/normal, plus a targeted panel sample with
`high_depth_mode`. Port the Java unit tests alongside each step, not at the end.

**6. Optimisation** - driven by profiles: single pass over the tumor BAM, threaded BGZF decoding,
flat arrays in place of hash lookups, phasing data structures. Parity suite runs on every change.

**7. Later** - TINC, append mode, then SBX and Ultima if needed.

## 6. Decisions and open questions

- Baseline sample: a WGS tumor/normal pair (REDUX BAMs). Paths live in `bench/env.local.sh`, which is
  not committed; outputs go to a directory outside the repository.
- The port lives in this fork under `sage-rs/`.
- Test data for SBX and Ultima: the public GIAB HG008 tumor/normal pair
  (`https://giab.s3.amazonaws.com/data_somatic/HG008/Liss_lab/`), sliced by chromosome and run through
  Java REDUX. REDUX 2.1 at this commit crashes during BQR because `BaseQualityResults.addPerfCounter`
  merges a counter into itself; the baseline was made with that one line patched locally. With
  `-form_consensus` REDUX also fails on the Ultima slices (a consensus read with an empty base), so the
  Ultima baseline was made without it.
