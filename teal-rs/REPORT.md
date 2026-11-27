# TEAL in Rust: record of the port

A record of how `teal-rs` was checked against Java TEAL 1.4 (this repository's `teal` module, Kotlin;
oncoanalyser runs TEAL 1.4 too), the faults found in the Java code, and timings. Options and outputs are
in [README.md](README.md).

## 1. How it was checked

Java TEAL 1.4 and the port were run on the same inputs with the same options, one Java run at a time.
- Java ran with `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2`. No TEAL output depends on
  identity hashes; the flag only fixes them.
- The port ran in Java mode (`TEAL_JAVA_COMPAT=1`).
- Both used 2 threads unless stated, which gives Java's prep one reader thread (TE6).

"Identical" means every output file is the same:
- **Telbam BAMs:** the header and every record, in order, compared as SAM text and as uncompressed
  BAM bytes (`scripts/bam_raw_equal.py`), bin, attribute order and attribute types included. Only the
  BGZF compression differs.
- **tellength.tsv:** byte for byte.
- **breakend.tsv.gz:** after decompression.

### Inputs

- a small private tumor/normal pair (chr21, chr22 and parts of chr1-10 of the first whole-genome pair;
  the tumor BAM is 2.8 GB, 65 million records);
- the public GIAB HG008 chr21 slices after REDUX: SBX (single-end reads with consensus records),
  Ultima, and the SBX CRAMs;
- the telbams of an earlier oncoanalyser run on a second private whole-genome pair (another person),
  written by Java TEAL 1.4, with that run's PURPLE, COBALT and BamTools files;
- synthetic BAMs (`scripts/make_synthetic.sh` and `make_synthetic.py`), described below;
- the port's own whole-genome telbams of the first whole-genome pair, with COBALT's Java output and
  BamTools' summaries for it, and the second pair's PURPLE purity. Only the depths, the duplicate
  proportion and the purity enter, and both tools got the same files, so the mix serves the comparison.

The synthetic BAMs contain:
- background pairs, and telomeric pairs with their mates before and after;
- split reads at shared and tumor-only break-end sites within a few bases of each other;
- SA supplementary reads;
- unplaced unmapped pairs in both mate orders, and unmapped placed mates;
- secondary and consensus (CR) records, duplicates and N bases;
- single-end reads and an alt contig;
- split reads at the GRCh37 blacklist boundary.

The set syn3 adds SA tags naming a contig that is not in the header.

### Comparisons

Times are wall, user CPU and peak RSS. The machine was shared with other jobs throughout.

| Case | Java | Rust | Result |
|---|---|---|---|
| Prep, small pair, tumor/normal | 218.1 s, 216 s, 1.57 GB | 17.3 s, 33 s, 0.21 GB | identical (2 files) |
| Prep, tumor only | 98.9 s, 100 s, 1.52 GB | 8.0 s, 15 s, 0.15 GB | identical |
| Prep, germline only | 110.8 s, 111 s, 1.53 GB | 9.0 s, 17 s, 0.14 GB | identical |
| Prep, tumor/normal, `-threads 1` | 195.7 s, 197 s, 1.57 GB | 33.9 s, 33 s, 0.12 GB | identical (and the same as at 2 threads) |
| Prep, HG008 SBX, `-sequencing_type SBX` | 86.4 s, 87 s, 1.54 GB | 8.5 s, 14 s, 0.15 GB | identical |
| Prep, HG008 Ultima, `ULTIMA` | 192.6 s, 190 s, 1.56 GB | 23.6 s, 39 s, 0.12 GB | identical |
| Prep, HG008 SBX CRAM, `-ref_genome` | 206.4 s, 212 s, 2.15 GB | 41.8 s, 46 s, 0.46 GB | identical (bytes too, with the CRAM attribute handling of section 2) |
| TealPipelineApp, the earlier run's telbams, tumor/normal | 133.3 s, 153 s, 2.73 GB | 11.6 s, 12 s, 0.71 GB | identical (5 files); also equal to that run's own tellength and breakend files |
| The same, `-threads 1` | 132.3 s, 150 s, 3.62 GB | 11.9 s, 12 s, 0.64 GB | identical; Java wrote 250 unplaced records twice (TE1) |
| The same, tumor only | 10.3 s, 21 s, 2.88 GB | 3.0 s, 3 s, 0.70 GB | identical |
| The same, germline only | 3.6 s, 11 s, 0.95 GB | 0.9 s, 1 s, 0.23 GB | identical; Java wrote 76 records twice (TE1) |
| The same, germline only, no `-purple_dir` | exit 1 | exit 1 | same failure (TE10); corrected, it runs |
| The same, `-ref_genome_version 37` | 121.9 s, 142 s, 4.43 GB | 11.0 s, 11 s, 0.71 GB | identical; Java: 88 records twice (TE1) |
| The same, `-sequencing_type SBX` | 132.0 s, 150 s, 4.35 GB | 11.5 s, 12 s, 0.70 GB | identical; Java: 260 records twice (TE1) |
| The same, `ULTIMA` | 131.3 s, 150 s, 3.84 GB | 11.7 s, 12 s, 0.70 GB | identical; Java: 114 records twice (TE1) |
| TealPipelineApp, the small pair's telbams | 3.9 s, 10 s, 1.07 GB | 1.0 s, 1 s, 0.08 GB | identical (5 files) |
| TealApplication (prep and analysis), depths on the command line | 211.6 s, 215 s, 1.63 GB | 20.1 s, 33 s, 0.20 GB | identical (5 files) |
| Synthetic sets 1 and 2: prep, TealPipelineApp, with GRCh37 and SBX (8 cases) | 1.0-1.4 s, 0.23-0.35 GB | 0.1 s, 0.03 GB | identical, all 8 |
| Synthetic set 3, SA tags naming an absent contig | exit 1 | exit 1 | same failure (TE11), same empty telbam left |
| TealPipelineApp, COBALT file without GC bucket 50 | 1.2 s, 0.34 GB | 0.1 s, 0.03 GB | identical (TE9 in Java mode) |
| TealApplication without GC50 depths | exit 1 | exit 1 | same failure (TE8) |
| TealPipelineApp, whole-genome telbams | 599.2 s, 627 s, 3.49 GB | 19.3 s, 19 s, 0.79 GB | identical (5 files); Java: 260 records twice (TE1) |
| Prep, first whole-genome pair, tumor/normal (96 and 120 GB BAMs) | 7,561 s, 7,271 s, 2.43 GB | 649 s, 1,232 s, 1.07 GB | identical (475,471 tumor and 331,014 normal records, uncompressed bytes); no duplicates from Java this time |

**TE1 in practice:** Java's duplicates are exact copies of unplaced unmapped records of complete groups,
adjacent in the sorted telbam. With them removed, Java's telbam equals the port's. They came in 6 of 9
Java runs on telbam inputs, and in none of the prep runs on full BAMs. No tellength or break-end output
changes, because TelbamReader drops a record its group already contains.

**Command line:** `-version` and `-help` of the three classes give Java's output byte for byte. 13
error cases give Java's exit status (1): missing ids or files, an unknown option, a bad version, thread
count or sequencing type, missing metrics and a bad depth. `teal -Xmx... <class> ...` works as
oncoanalyser calls it, and `teal -version | sed -n '/Teal version/ { s/^.* //p }'` gives 1.4.

**Corrected mode:** the telbams never change, as no switch touches the prep except TE11. On the
whole-genome telbams:
- telomere length moves -0.08% (normal) and -0.04% (tumor, raw and final), as about 0.1% of the
  fragments change class (TE2);
- of 21,059 break ends, 220 rows change: about 70 move position and 82 lose the "duplicate" filter
  (TE4);
- the 84 PASS break ends are the same.

**Unit tests:** `cargo test --release` runs 29 tests:
- the ported Kotlin tests: TealUtilsTest, ReadGroupTest, TelomereMatcherTest, TelomereLengthCalcTest,
  BreakEndSupportCounterTest and TelomericSplitReadAnalyserTest;
- the port's own tests:
  - the packed telomere test against the text test;
  - the aligner's tie and gap rules;
  - FragmentClassifier's regex backtracking;
  - read group acceptance and missing regions;
  - region consolidation and order;
  - htsjdk's integer types;
  - option parsing;
  - an end-to-end prep on a BAM written by the test, at 1, 2 and 3 threads.

**Whole-genome prep:** on the first whole-genome pair at `-threads 2`, both telbams are identical to Java's,
uncompressed byte for byte. That is 475,471 tumor and 331,014 normal records. The port took 649 s against
Java's 7,561 s, 11.6 times faster, and 1.07 GB against 2.43 GB.

### In the pipeline (GIAB HG008, 2026-11-27)

The port in Java mode ran TEAL_PREP and TEAL_PIPELINE on Google Batch, against the Java TEAL 1.4 of the
first HG008 run on the same inputs.
- **Lengths and break ends:** both tellength files and the breakend file are identical.
- **Telbams:** the same records in the same order (165,022 tumor, 104,103 normal), although Java's prep ran
  at 16 threads (TE6). The headers differed only in `@HD VN`: Java writes htsjdk's version, 1.6, and the
  port kept the input's 1.5. From 2026-11-27 the port writes 1.6 too (hmf-common-rs `sam_header`). A
  `VN:1.5` slice of the HG008 tumor BAM then gives a byte-identical telbam header from both. Every local
  input said 1.6.
- **Unplaced records:** after that fix, the telbams were identical as SAM text but not as BAM bytes.
  11,792 unplaced records had the input's bin, 4680, where htsjdk writes 0 for any record without a
  reference. The port now writes 0 for them too. On a public HG008 test BAM (the start of chr1 and
  200,000 unplaced reads), Java and the port give a byte-identical telbam (uncompressed, 1,916 records).
  HG008 rerun 7 (2026-11-27) confirmed it in the pipeline: both whole-genome telbams are byte-identical
  to the first run's Java telbams once decompressed (66.3 and 42.1 MB). The prep took 6 min 52 s at 439%
  CPU and 1.5 GB.
- **Reading on Batch:** the read-ahead (several handles, in file order) took the prep from 15.5 to 7.5
  minutes at 420% CPU, against Java's 13.8 (HG008 rerun, 2026-11-27).
- **Time:** TealPipelineApp took 7.4 s and 0.4 GB against Java's 35 s and 2.3 GB. The prep took 15.5
  minutes at 165% CPU and 1.1 GB, against Java's 13.8 minutes at 1112% CPU and 8.6 GB (n2-highmem-16,
  BAMs through gcsfuse). The port reads the two BAMs one after the other, each as one stream, and one
  gcsfuse stream gives about 160 MB/s. Java's 15 reader threads read in parallel (section 5).

## 2. How the prep reproduces Java

BamProcessor runs `max(threads - 1, 1)` reader threads. With one, at `-threads` up to 2, the writer sees:
1. the unplaced unmapped records;
2. every placed record in file order (10 Mb partitions of each contig, in header order);
3. in each missing-read pass, the unplaced records again and the consolidated mate and supplementary
   regions, until a pass accepts no read. The regions are sorted by chromosome rank, then name, then
   start, and merged within 100 kb.

The port reads in exactly that order at every thread count, and runs the same state machine:
- a telomeric record (TTAGGGTTAGGG or CCCTAACCCTAA; CR consensus records skipped) opens a group;
- records with a group's name join it;
- a group with both reads and all its SA-listed supplementary reads is written and closed;
- secondary alignments join as supplementary, with their flags changed as Java's acceptRead does;
- no group opens in the missing-read passes;
- the incomplete groups are written at the end in Java HashMap order.

The records are then sorted as htsjdk's writer sorts them (SAMRecordCoordinateComparator, stable), so
the order they were added in only decides ties. The header is the input's.

**CRAM records as htsjdk gives them.** The port does the same as htsjdk to the records it keeps from a
CRAM, in both modes, since this is library behaviour rather than a fault:
- MD and NM are not regenerated;
- attributes are kept sorted by the tag's binary code;
- integers are written with the smallest type;
- bin 0 is written for records without a position.

From a BAM, both write the records' bytes unchanged.

**The analysis.** TealPipelineApp runs the prep again on the telbams it is given and writes new
telbams; the lengths and break ends come from those. On the earlier run's telbams, the new telbams equal
the inputs, in Java and in Rust. The break ends sit in an arena, so that Java's sharing of
TelomericBreakEnd objects between the tumor and germline analyses (TE4) can be reproduced.

## 3. Findings

Corrected by default, with Java's behaviour from `-java_compat` or `-java_compat_except TE<n>,...`:
TE2, TE4, TE7, TE8, TE9, TE10 and TE11. No switch:
- TE1 and TE6 depend on timing, and the port gives Java's single-reader, race-free result in both
  modes;
- TE3 and TE5 keep Java's behaviour in both modes, as the intent is unclear.

| Id | Java | Effect | Port |
|---|---|---|---|
| TE1 | A race: BamProcessor submits the tasks of a missing-read pass before it calls `writer.setProcessingMissingReadRegions(true)`. The first task, the unplaced reads, can start before the flag is set. A telomeric unplaced record whose group was completed in the first pass then opens a new group, which takes its mate again, so the fragment is written twice. | 76 to 260 unplaced records written twice in 6 of 9 Java TealPipelineApp runs; the counts differ for the same input (0 and 250 at 2 and 1 threads). Lengths and break ends are unchanged. | The flag is set before the pass, in both modes. |
| TE2 | TealUtils.reverseComplementSequence appends nothing for N or any other non-ACGT base, so the result is shorter and shifted: an N inside a telomeric run can join two partial hexamers or split one. Used by FragmentClassifier and BreakEndSupportCounter. | On the earlier run's telbams, 0.1-0.3% of the F1/F2/F4 counts change; telomere length moves -0.11% (tumor) and -0.14% (normal). Break ends are unchanged. | N stays N, and IUPAC codes are complemented. |
| TE3 | TelomereMatcher computes its match ratio (with integer division, so 0 or 1) and only logs it: any alignment of 12 matched bases counts. The 0.9 and 0.8 thresholds of BreakEndParams have no effect, unlike the README's "at least 90% match ... at least 80% match". | | Kept: which ratio was meant (over the aligned length or the clip) is unclear, and applying one changes the candidates. |
| TE4 | The germline analysis starts from the tumor's TelomericBreakEnd objects, and the data class `copy()` keeps them shared. Merging candidates within 2 bases then moves the tumor's break end after its support was counted, and markDuplicateBreakEnds marks the shared object duplicate for both. | On the earlier run's telbams, 87 of 4,965 break ends move to another position (the tumor's own), 87 lose the "duplicate" filter and 2 gain it. On the whole genome, about 70 move and 82 lose the filter; the 84 PASS calls are unchanged. | The germline analysis gets copies. |
| TE5 | populateLongestTelomereSegment searches the mate as stored, while classifyPairedRead first turns it to point away from the aligned read. A mate classified "not telomeric" can still give a 150-base maxTelomericLength, and the reverse. | With the mate oriented, 1,690 of 4,965 break-end rows would change (PASS unchanged). | Kept: the README defines maxTelomericLength without an orientation. A question for Hartwig. |
| TE6 | With `-threads` 3 or more, several readers feed one writer, so which records reach a group before it completes depends on thread timing. | Output depends on timing. | The single-reader result at every thread count. |
| TE7 | loadExcludedRegionBed copies the GRCh37 blacklist to a temporary file whose name does not end in `.bed`, so its 0-based starts are read as 1-based. | GRCh37 only: a split read at a region's BED start lies inside the region, so its candidate is dropped, as is a read group lying entirely within such regions. Synthetic GRCh37 case: one more break end when corrected. | Start + 1. |
| TE8 | `-reference_gc50_read_depth` and `-tumor_gc50_read_depth` are optional, but TealApplication passes null to calcTelomereLength (`gc50ReadDepth!!`: NullPointerException, exit 1). | The run fails. | Falls back to the mean depth, as TelLengthApp.run does. |
| TE9 | CobaltGcMedianFile gives -1 for a bucket it does not list, and TealPipelineApp takes that as the GC 50 read depth: the GC bias becomes -1/mean and the adjustment 1.2427. | Synthetic case: lengths 29% higher than with the mean. | The mean depth, as in TE8. |
| TE10 | TealPipelineApp registers `-purple_dir` as required, although validate() needs it only with a tumor. | A germline-only run without it fails. oncoanalyser passes `-purple_dir` only when it has a PURPLE directory. | Optional. |
| TE11 | A missing mate or SA position on a contig that is not in the header is queried with htsjdk, which throws "Invalid reference index -1". | The whole run fails (exit 1, an empty telbam). This happens with SA tags written against another reference, or with a BAM re-headered with fewer contigs. | Such positions are skipped with a warning. |

Kept as Java, being either not faults or of unclear intent:
- TealPipelineApp makes the telbams again from the telbams it is given, so the prep's output is read
  twice.
- The README differs from the code in places: duplicates are called within 60 bases, not "50 bases"
  (`markDuplicateDistance`), and discordant pairs within 500 bases, not "1000 bases".
- A record whose CIGAR covers no reference base, starting at a partition start, is missed by htsjdk's
  overlap query, and so by Java. Such records are invalid for mapped reads and did not occur here.
- Long CIGARs stored in a CG tag are not expanded by the port's reader; none occur in short-read data.

## 4. Timings and memory

**Small pair and slices** (section 1): at 2 threads Rust is 8 to 13 times faster than Java on BAM and 5
times faster on CRAM, which htslib decodes. It uses 5 to 13 times less memory.

**Scaling:** prep of the small tumor, with byte-identical output each time.

| -threads | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| wall | 17.7 s | 7.9 s | 5.3 s | 4.0 s |
| user CPU | 15.3 s | 15.1 s | 15.0 s | 15.0 s |
| peak RSS | 0.08 GB | 0.19 GB | 0.21 GB | 0.29 GB |

Java takes 98.9 s and 1.52 GB at 2 threads (one reader thread). The port's CPU time does not grow with
threads: every thread inflates BGZF batches with libdeflate, so the work is shared, not duplicated.

**Whole-genome prep** of the first whole-genome pair: tumor 96 GB, normal 120 GB, with 54.5 and 51.2 million
unplaced reads; `-threads 2`, Java mode. With the read-ahead (4 handles, 32 MB chunks) it took 649 s wall
for both BAMs, 1,232 + 72 s CPU, and 1.07 GB peak; before the read-ahead, 827 s and 0.83 GB. The extra
memory is the chunks read ahead. Java took 7,561 s, 7,271 + 252 s CPU and 2.43 GB on the same machine,
with the same telbams. In an oncoanalyser run on a pair of that size, the Java prep used about 5 CPU
hours and peaked at 18.5 GB.

Where the time goes, from macOS `sample` on the normal's main pass:
- inflating (libdeflate, with its CRC32): about 60%;
- the telomere test: 11%;
- the kernel's read copies: 11%;
- the first touch of each record: 7%.

For comparison, `samtools view -c` of the small tumor takes 18.4 s with `-@1` and 6.1 s with `-@3`,
against the port's prep at 7.9 s (2 threads) and 4.0 s (4 threads).

**TealPipelineApp** takes 11-12 s against Java's 122-133 s on the earlier run's telbams (0.7 GB against
2.7-4.4 GB), and 19 s against 599 s on the whole-genome telbams. The support counter looks only at
groups with a read within 500 bases of a break end, and caches each read part's telomere alignment.
Java tries every group against every break end.

**Memory:** the prep holds the open groups and the telbam's records until they are sorted (0.3-0.5 GB
at whole-genome sizes), plus up to 2 x threads + 2 inflated batches of about 4 MB each.

## 5. Left out, and open points

- **Stand-alone programs not ported:** TelbamApp, TelLengthApp and BreakEndApp. Their code is ported
  as far as the pipeline apps use it, but not their command lines, so there is no `-tsv_file`,
  `-specific_chr` or `-genome_regions`.
- **TE3 and TE5** are questions for Hartwig.
- **Not run:** Java with 3 or more threads (TE6 says why its result varies), and a whole-genome CRAM.
- **Prep reading on Batch:** the read-ahead (several ranged reads ahead of the record the reader is at)
  took the prep from 15.5 to 7.5 minutes. Reading the tumor and normal at the same time could take it
  further. It is not done, because 7.5 minutes is already under Java's 13.8.
