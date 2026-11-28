# markdups-rs: how it was checked, Picard's faults, timings

Picard's MarkDuplicates as GATK 4.6.1.0 runs it (Picard 3.3.0, htsjdk 4.1.3). Made 2026-11-28.

## Checks (all under -java_compat; `scripts/regress.sh`, log in `baseline/markdups/regress.log`)

Comparison (`compare.sh`): `samtools view --no-PG -h` of both BAMs, the uncompressed BAM bytes (`gzip -dc`),
the metrics file without its "Started on" line, the exit code. Every case runs from a folder holding `<id>.bam`
so that GATK and the port record the same command line. GATK references: the pipeline's image through docker
(`--platform linux/amd64 --cpus N`, `-Xmx4000M`), and the same jar on the Mac's Java 23 for synthetic cases; on
the RNA slice docker and jar outputs are identical (records and uncompressed bytes; the jar uses the JDK
deflater, the image Intel's: compressed bytes differ, nothing else).

| Case | Records | Duplicates | GATK reference | markdups-rs -java_compat, 1 thread | 2 threads |
|---|---|---|---|---|---|
| RNA slice (TUMOR) | 3,278,190 | 1,956,117 | docker 1 and 2 CPUs; jar 1 and 2 threads (all four identical) | identical | identical |
| HG008-T HLA slice | 2,197,409 | 368,882 | docker 1 and 2 CPUs | identical | identical |
| HG008-N-D HLA slice | 1,216,234 | 142,288 | docker 1 and 2 CPUs | identical | identical |
| HG008-T IGH slice | 1,010,206 | 170,493 | docker 1 and 2 CPUs | identical | identical |
| One step, RNA slice: unsorted pair-grouped input | 3,278,190 | 1,956,117 | samtools 1.21 sort then GATK, both docker | identical (in memory, and with 10 runs on disk) | identical (both) |
| 30 synthetic cases | 0 to 140,000 | | jar | identical | identical |
| 6 synthetic cases GATK refuses (unknown read group, invalid record, unsorted records, unsorted header, no @HD, MD4) | | | jar | refused with the same message | same |
| 3 cases of several inputs (read group renamed, identical headers, three files) | 1,200 to 1,800 | | jar | identical | identical |

`regress.sh` (2026-11-28): 0 differences; 78 comparisons identical, 12 identical refusals. "Identical" means the
decoded records and header (`samtools view --no-PG -h`), the uncompressed BAM bytes and the metrics (but the date
line) all equal.

Index: the port's `.bai` against htsjdk's `BuildBamIndex` of the port's own BAM (`bai_diff.py`): identical but for
the last record's chunk end and the last reference's metadata end offset, which htsjdk's writer (CREATE_INDEX)
puts inside the last data block and BuildBamIndex at the start of the end-of-file block; GATK's own `.bai` for
its own BAM differs from BuildBamIndex in the same two places (unit test `bai::tests::as_htsjdk_writer` holds
GATK's bytes). samtools reads the index (idxstats, region queries).

Inputs:
- RNA slice (private, sample id TUMOR): `baseline/isofox/inputs/TUMOR.rna.slice.bam` is MarkDuplicates output;
  `make_rna_input.sh` removes the `PG:Z` tags and 0x400 flags and gives it the header samtools sort writes after
  STAR (`@HD VN:1.4 SO:coordinate`, STAR's and samtools' @PG lines with generic command lines). The record order
  is samtools sort's (checked: (tid, pos, reverse) never decreases). 3,278,190 records; 3,281 pairs have their
  mate outside the slice (they stay unmatched, as in Picard).
- Unsorted, pair-grouped version for the one-step mode (`make_unsorted.py`): the RNA input sorted by name
  (samtools sort -n), each read's records in STAR's Unsorted order (alignment by alignment in HI order, the
  leftmost mate first, supplementary records after their alignment), header as STAR writes it (`@HD VN:1.4`, no
  SO, no samtools @PG). Close to STAR's real output in grouping and within-read order; the order of reads differs
  (STAR writes them in FASTQ order, here name order) and the tags are in htsjdk's order (STAR's own order is
  covered by the synthetic `tags` case). The comparison does not depend on it: the reference is samtools 1.21
  sort of this file followed by GATK.
- HG008 (public): the tumour and normal HLA slices and the tumour IGH slice (REDUX output; NovaSeq names, 8 read
  groups in one library, supplementaries, unmapped reads): REDUX's consensus records (CR tag) removed, 0x400
  cleared, the slicing's @PG lines (with storage paths) dropped (`make_hg008_inputs.sh`).
- Synthetic (`make_synthetic.py`, 36 cases, plus three of several inputs): pairs in all orientations, soft and hard clips, spliced reads,
  fragments, mates on other references, secondary and supplementary records, unmapped mates and unplaced pairs,
  four read groups in two libraries plus "Unknown Library", optical sets of 2 to 30 with 5- and 7-field names,
  unparseable names and keepers without location, score ties, missing qualities, QC failures, empty input,
  unsorted records, unsorted and missing @HD, STAR's tag order and types (arrays, floats, chars, values above
  2^31), existing DT and PG tags, a 70,000-pair set (Java's HashMap order of the optical cluster map), every
  option the port implements, invalid records under STRICT and LENIENT, the faults below, and several inputs (RG
  collision renamed, identical headers, three files).

Unit tests: `cargo test --release`: 51 passed (markdups-rs); `hmf-common-rs`: 143 passed, 1 ignored. Picard's tests of the ported classes are in `src/picard_tests.rs`
(the AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesTest cases, through a tester built like
Picard's SamFileTester and MarkDuplicatesTester), `src/optical.rs` (OpticalDuplicateFinderTest,
ReadNameParserTest), `src/metrics.rs` (library sizes against Picard's own `DuplicationMetrics.main`).

## Findings (MD1 to MD7)

Each is corrected by default; `-java_compat` (or `-java_compat_except`) reproduces Picard 3.3.0.

**MD1: optical-duplicate coordinates kept as 16-bit numbers.** ReadEnds extends PhysicalLocationShort, whose
setX/setY cast to `short`; ReadNameParser also casts the tile to `short`. Cluster coordinates above 32767 wrap
(32768 becomes -32768). Wrong because reads 1 pixel apart across 32767/32768 are 65535 apart for the distance
test, and coordinates 65536 apart look identical; and because Picard orders tied duplicates by tile, x, y, the
wrapped values also change which read of a tied set is kept. FastDup reports the overflow and keeps it on purpose.
Evidence: the synthetic `md1_coords` case (identical to GATK under -java_compat; corrected, READ_PAIR_OPTICAL_
DUPLICATES 4 -> 3 and ESTIMATED_LIBRARY_SIZE 6 -> 4, the pairs judged on their real coordinates). Effect on the RNA slice (patterned flowcell;
12% of the reads have y > 32767): corrected, 3 tied duplicate pairs keep a different pair (12 records change
flag, 6 gain and 6 lose 0x400), READ_PAIR_OPTICAL_DUPLICATES 19,851 -> 19,885, ESTIMATED_LIBRARY_SIZE
558,144 -> 558,147, the optical/non-optical set histograms change slightly. Corrected: full int coordinates.

**MD2: read-name numbers that do not fit an int wrap.** ReadNameParser.rapidParseInt accumulates `val * 10 + d`
in an int without a check, so a field of more than 10 digits becomes an arbitrary coordinate, which can be
"close" to an unrelated read's. Corrected: such a name has no location (as a non-numeric field). Evidence:
`md2_numbers` (GATK counts 2 optical duplicates among reads whose 12-digit fields wrap close together;
corrected 0). Effect: none on Illumina names.

**MD3: MAX_OPTICAL_DUPLICATE_SET_SIZE -1 turns optical detection off.** The option's documentation: "To
completely disable this check, set the value to -1". OpticalDuplicateFinder's constructor stores the value as
given (only the unused setter maps values below 1 to Long.MAX_VALUE, and then overwrites it with the argument)
and findOpticalDuplicates returns no flags when `length > maxDuplicateSetSize`, which is every set for -1.
Corrected: values below 1 mean no limit. Evidence: `md3_maxset` (GATK: 0 optical duplicates with -1;
corrected: 3, as with any limit above the set size).

**MD4: paired reads without a known read group crash Picard.** buildSortedReadEndLists makes the mate key from
`rec.getReadGroup().getReadGroupId()`; for a paired read with a mapped mate and no RG tag that is a
NullPointerException (GATK exits with it). Corrected: the key uses the RG value (empty without one). Evidence:
`md4_no_rg` (GATK: `java.lang.NullPointerException: Cannot invoke ... getReadGroupId()`; the port under
-java_compat stops with the same message).

**MD5: the mate key joins read group id and read name without a separator.** `key.append(rgId);
key.append(readName)`, so read group "A" with read "Bx" and read group "AB" with read "x" share a key and are
paired as mates; the pair ends, and so the duplicates, are then wrong. Corrected: a NUL byte between them.
Evidence: `md5_rg_name`: two pairs whose first mates GATK pairs with each other across read groups, and a real
pair at those two positions; GATK (and the port under -java_compat) marks the real pair as a duplicate of the
made-up one, corrected it is not marked. Effect: needs read group ids that are prefixes of others with names that
end like others; none in Illumina data.

**MD6: without a read name regex the set-size histograms count only single pairs.** trackDuplicateCounts (which
fills all_sets and non_optical_sets) is called from trackOpticalDuplicates, which MarkDuplicates skips when
READ_NAME_REGEX is null, while addSingletonToCount runs for every chunk of one pair. The histograms then claim that
every duplicate set has size 1. Corrected: every chunk counts, with no optical duplicates. Evidence:
`opt_no_regex` (its two sets of 5 pairs are missing from GATK's histogram). Effect: metrics histograms only, and only with READ_NAME_REGEX null.

**MD7: PG tags of the unplaced reads after the first can be removed.** Pass 1 stops at the first unmapped read
without a reference ("no reason to continue"), so PG values that occur only on later reads are never collected;
pass 2 sets `PG` to `chainedPgIds.get(oldValue)`, null for those, which removes the tag. Corrected: the values
of all records are collected (the port reads them anyway for its sort-order check), so they get their @PG lines
and keep a tag. Evidence: `md7_pg_tail` (corrected: 4 reads keep `PG:Z:MarkDuplicates.3`, chained by a second
@PG line to samtools). Effect: none for STAR output (no PG tags on records).

### Other observations (behaviour kept)

- Unmapped mates of mapped reads and secondary/supplementary records are never marked in coordinate mode
  (documented by Picard; ISOFOX relies on primary flags).
- Picard's results do not depend on its heap: the ReadEnds codecs used when spilling round-trip every field.
- htsjdk's chunk merging in the index tests `block1 + 1 == block2` on block addresses (bytes), so only chunks in
  the same block are merged; kept, for the same index.
- htsjdk's out-of-order error names the offending record twice ("Offending records are at [chr1:100] and
  [chr1:100]"); the port names both records. Both fail.
- Barclay prints a given null (`--READ_NAME_REGEX null`) as an empty word, so the @PG CL has two spaces; kept.
- ROI histogram values above 1 on slices (CoverageMult 1.0026...): READ_PAIRS_EXAMINED counts reads whose mate is
  outside the file; not a fault for whole BAMs.
- The JVM with `-Xmx` of 80% of the task's memory (the pipeline) peaks at most of it (59 GB seen); not Picard's
  need: the slices ran in under 3 GB with -Xmx4000M.

## Timings and peak memory (`/usr/bin/time -l`; docker: the container's cgroup `memory.peak`)

Measured on this Mac (arm64, shared with other delegates' runs; some of these ran beside another Java run). The
docker GATK runs execute the amd64 image under emulation, so the native jar (Java 23, same records) is the fairer
Java figure; GATK's heap was `-Xmx4000M` (the pipeline gives 80% of the task's memory and peaked at 59 GB).

| Case | GATK docker 1 CPU | GATK docker 2 | GATK jar 1 thread | GATK jar 2 | markdups-rs 1 thread | markdups-rs 2 |
|---|---|---|---|---|---|---|
| RNA slice | 39.4 s, 2655 MB | 32.2 s, 2993 MB | 18.8 s, 2500 MB | 20.6 s, 2583 MB | 7.5 s, 208 MB | 5.0 s, 297 MB |
| HG008-T HLA | 37.5 s, 2750 MB | 21.1 s, 2661 MB | | | 5.5 s, 313 MB | 3.7 s, 397 MB |
| HG008-N-D HLA | 21.8 s, 2432 MB | 14.7 s, 2386 MB | | | 3.0 s, 181 MB | 2.0 s, 233 MB |
| HG008-T IGH | 20.9 s, 2395 MB | 13.2 s, 2430 MB | | | 2.7 s, 126 MB | 1.8 s, 187 MB |

Against the native jar on the RNA slice: 2.5 times faster at 1 thread, 4.1 times at 2, a twelfth of the memory.
GATK gets nothing from a second thread (its MarkDuplicates is single-threaded; the second CPU goes to the JVM).
The port's output BAM is 8% smaller than GATK's at the same level (libdeflate vs Intel's deflater: 164 vs 178 MB).

One step, RNA slice (the unsorted input of 3.28 M records):

| | time | peak memory | disk |
|---|---|---|---|
| samtools 1.21 sort, docker, 1 thread | 19.8 s | 1329 MB | sorted BAM 134 MB (level 6) |
| samtools 1.21 sort, docker, 2 threads | 9.1 s | 1834 MB | same |
| then GATK MarkDuplicates, docker, 2 CPUs | 28.5 s | 3174 MB | |
| two steps, 2 threads, in all | 37.6 s | 3174 MB | 134 MB kept until the task ends |
| SortAndMarkDuplicates, 1 thread, SORT_MEMORY_MB 2048 (all in memory) | 8.6 s | 1480 MB | none |
| SortAndMarkDuplicates, 1 thread, SORT_MEMORY_MB 128 (10 runs) | 12.7 s | 344 MB | runs at level 1, about 1.3 times the sorted BAM, deleted at the end |
| SortAndMarkDuplicates, 2 threads, SORT_MEMORY_MB 2048 | 7.0 s | 1482 MB | none |
| SortAndMarkDuplicates, 2 threads, SORT_MEMORY_MB 128 | 10.4 s | 344 MB | as above |

Against the port's own two steps (samtools sort, then markdups-rs MarkDuplicates at 5.0 s), the one step saves the
sort's level-6 compression of the whole BAM, its writing, the two reads of it, and the disk it holds; the price is
the sort buffer (SORT_MEMORY_MB plus 32 bytes per held record) or level-1 runs read twice.

Profile (macOS `sample`, RNA slice, 1 thread, 8.6 s): deflating the output (level 2, as GATK) 42% of the
samples, inflating the input twice 10%, parsing and re-encoding tags 17%, sorting the ends 4%, the rest record
handling. More threads go to deflating (pass 2) and inflating (pass 1).

## Whole-sample estimate

The RNA slice is 3.28 M records (about 2% of the 74 M-pair tumour RNA, some 150 M records): about 4 min at 2
threads. More threads speed up inflating and deflating, but each pass handles the records on one thread (pass 1
at about 1.3 M records/s on one thread with inflating, pass 2 of the same order), so about 2 min at 16 threads; splitting that work by
reference is the next step if needed. GATK took 15.3 min there. Memory is one reference's ends (about 40 bytes
per read, 56 per pair) plus the pending mates and a bit per record: about 1-2 GB for the largest chromosome of that RNA. SortAndMarkDuplicates
holds SORT_MEMORY_MB of records plus 32 bytes per held record, writes runs at level 1 to TMP_DIR (about a
quarter of the uncompressed size), and saves the sorted BAM's level-6 compression, its write and two reads, and
its disk space.


## The whole sample in the pipeline (Google Batch, 2026-11-28)

The 74 M read-pair tumour RNA, through oncoanalyser's STAR route (fork profile `rust_ports`, `--rust_java_compat`;
STAR 2.7.3a and samtools sort from the cache), image `55ce55b212`, against GATK 4.6.1.0's task of the first RNA run on
the same sorted BAM:

| | GATK 4.6.1.0 | markdups-rs, `--threads 16` |
|---|---|---|
| machine | n2-standard-64 (the alignment selector; the pipeline's -Xmx is 80% of 200 GB) | n2-standard-16 |
| wall time of the tool | 15 min 19 s | 4 min 45 s (286% CPU) |
| peak memory | 59.4 GB | 1.7 GB |
| VM cost | $0.98 | $0.09 |
| BAM size | 9.08 GB | 8.38 GB |

- Records: `samtools view --no-PG -h` of the two BAMs has the same MD5 (164.7 M records, 99.3 M marked duplicates).
- Metrics: identical apart from the `# Started on` line.
- Index: `samtools idxstats` through each BAM's own `.bai` identical, and region counts on four regions equal.
- Downstream, against the earlier rust_ports run on GATK's BAM: ISOFOX 2.0.1 mode (all 9 files), CIDER, CUPPA and
  BamTools byte-identical, apart from BamTools' `partition_stats` (row order and process times, as between any two runs).

The time is mostly the two passes' single-threaded record handling (the estimate above); the memory is one
reference's ends, as estimated.
