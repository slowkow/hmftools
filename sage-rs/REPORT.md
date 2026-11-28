# sage-rs: progress record

A running record of what was built, how it was checked, how it performs and what was found in the
Java code along the way. It is written so that it can later be turned into an announcement, an
upstream issue list or a pull request description. Detail lives in [PLAN.md](PLAN.md) (status and
measurements) and [QUIRKS.md](QUIRKS.md) (the Java findings by subject).

Baseline: Java SAGE 5.1 at upstream master `883ab11bba` (2026-11-27); until then `a224a2c` (`sage-v5.1-beta.4`).
Section 9 has the update from one to the other.

Two modes. By default the port corrects 33 of the 60 Java findings of section 5. With `-java_compat`
(or `SAGE_JAVA_COMPAT=1`) it behaves as Java does. **Sections 1, 2, 3 and 6 describe the port run with
`-java_compat`**: every "identical", every record count and every timing there was measured in that
mode, or before the corrections existed, when it was the only mode. Section 5.5 has what the corrected
mode changes, and its time and memory against Java mode.

## 1. What exists

`sage-rs` is a Rust reimplementation of SAGE 5.1. It takes the same command-line options and, with
`-java_compat`, writes the same VCF, byte for byte, as Java run with `-threads 1`.

| Area | State |
|---|---|
| Tumor/normal, tumor-only, several tumors / references | identical to Java |
| `-panel_only`, `-high_depth_mode`, `-germline`, `-disable_soft_filter` | identical |
| Illumina, SBX, Ultima | identical |
| MSI jitter parameters, split (RNA) reads | identical on simulated MSI input; split reads have unit tests only |
| TINC (`-run_tinc`, standalone `sage-rs tinc`) | identical, including simulated contamination at 2%, 10%, 50% |
| Append (`sage-rs append`) | identical |
| CRAM 3.0 | identical to Java reading the same CRAM |
| CRAM 3.1 | works; Java cannot read it |
| Visualiser (`-vis_*`) | same files as Java once the read rows Java orders at random are sorted (section 6) |
| Unit tests | the Java SAGE unit tests are ported, plus tests for the port's own code |

Section 6 lists what was added after this table was first written.

## 2. How it was checked

- **Definition.** The decompressed VCF of the port run with `-java_compat` must equal Java's with
  `diff`. Java is run with one thread because its output depends on thread scheduling (section 5.2);
  the port gives the same bytes at any thread count, in both modes (in Java mode with the exception
  of S24, where it follows Java).
- **Data.** A private Illumina WGS tumor/normal pair (whole genome), and the public GIAB HG008
  tumor/normal pair for SBX and Ultima (chr21, after Java REDUX).
- **Whole genome.** 360,440 records without TINC and 611,876 / 365,286 with `-run_tinc`. Against Java
  at 16 threads the only differing lines are those with thread-order `LPS` ids (4,573 and 4,704).
- **A second private project** (run on the pipeline's own REDUX and PURPLE output; Java 5.1 run fresh,
  since the pipeline had used SAGE 4.2, with 4.2's recalibration table under the 5.1 header):

  | Input | Run | Records | Result |
  |---|---|---|---|
  | Tumor RNA of three samples, appended to the PURPLE somatic and germline VCFs (`-skip_bqr -skip_msi_jitter -max_read_depth 100000 -write_frag_lengths`) | Java, 1 thread | 4,248 to 11,021 somatic, 1,910 to 1,957 germline | all six VCFs and fragment length files byte-identical |
  | Deep exome pair, whole run | Java, 16 threads | 1,198,874 | 97 lines differ, `LPS` ids only |
  | WGS pair from another vendor | Java, 16 threads; chr21 at 1 thread | 430,927; 3,751 | 429 lines differ, `LPS` ids only; chr21 byte-identical |
  | WGS pair whose tumor has no base quality above 28 | same | 1,945,878; 13,923 | every line identical; chr21 byte-identical |
  | WGS pair, a later sample | same | 390,056; 4,707 | 2,094 lines differ, `LPS` ids only; chr21 byte-identical |

  Read from a network drive (about 260 to 400 MB/s), the three whole genomes took Java 3,531, 4,132 and
  4,128 s with 18 to 22 GB, and the port 620, 1,226 and 862 s with 1.5 to 1.8 GB; the exome 1,549 s
  against 383 s. The port was waiting on the drive for most of that time.
- **Regression suite.** 49 comparisons against stored Java output, run in Java mode after every
  change: each mode on each platform, MSI, multi-sample, append, CRAM, fragment lengths, specific
  variants, TINC with two tumors. The same inputs are then run in the corrected mode (section 5.5). It is kept outside the repository because the
  outputs contain patient variants; `bench/` has the scripts that produce each side.
- **Unit tests.** One Rust module per Java test class in `src/tests`.

## 3. Performance

Same machine (20 cores), BAMs on an external volume. One thread, chr21:

| Input | Java | Rust | Ratio |
|---|---|---|---|
| Illumina BAM | 501 s, 3.6 GB | 59 s, 150 MB | 8.5x |
| SBX BAM | 106 s | 37 s | 2.9x |
| Ultima BAM | 253 s | 109 s | 2.3x |
| Illumina CRAM | 767 s | 110 s | 7.0x |
| SBX CRAM | 310 s | 81 s | 3.8x |
| Ultima CRAM | 692 s | 201 s | 3.4x |

Whole genome, 16 threads:

| Run | Java | Rust |
|---|---|---|
| Without TINC | 2,450 s, 26.9 GB | 812 s, 0.9 GB (first complete version) |
| With `-run_tinc` | 3,407 s | 764 s (first complete version), 593 s (second round), 327 s (third), 316 s and 1.7 GB (current) |
| Standalone TINC | 113 s, 17 GB | 14 s, 0.5 GB |

The current whole-genome run uses 4,030 core-seconds (9,030 for the first complete version) and
averages 12 of 16 cores. Java's 3,407 s was measured while other jobs were running, so the ratio
(about 10x) is approximate.

### What made it faster, in order of effect

1. A straight port: no per-mismatch `String`s, no `SAMRecord` built per merged fragment, array lookups
   in place of string-keyed hash maps.
2. The tumor's reads are read once per partition and kept for the evidence stage (Java queries the
   file twice).
3. BAM tags parsed in place; read buffers reused.
4. Overlapping mates that are both a plain match are merged without the cigar walk.
5. The candidate stage skips from mismatch to mismatch with a word-wise compare.
6. A 1 MB read buffer per alignment file handle. With htslib's default block the 16 threads made many
   small requests to the disk and the whole genome read at 340 MB/s; with 1 MB blocks it takes 327 s
   in place of 593 s. This is the largest single gain in wall time and changes no code path in the caller.

### Tried and dropped

- **mimalloc as the default**: 16% less CPU time, but 11.6 GB against 1.7 GB on a whole genome, and
  only 5% less wall time there (301 s against 316 s) because the external disk is the limit. Its purge
  settings and a collect on each worker between regions did not lower the memory. It is kept as a
  build feature (`--features mimalloc`), off by default.
- **Slice comparison in the repeat search**: slower than the byte loop for the short repeats involved.

### Where the time goes now

Illumina chr21, one thread: BAM decompression and reading 34%, matching reads to candidates and
quality scores 13%, allocation 11%, decoding records 8%, jitter check 7%, read-context building 6%,
fragment sync 5%. For Ultima, read-context building is about half.

At 16 threads the whole genome is limited by the disk, not the CPU. From an internal SSD (chr1 and
chr2) the run scales with the threads (82 s on 8, 43 s on 16) and keeps 15.4 of 16 cores busy; from the
external disk the whole genome reads 690 MB/s of the 800 the disk can give, so about 270 s is the floor
there. The binary already uses libdeflate. Busy time at 16 threads: decompression 22%, allocation 10%,
copying 8%, record decoding 7%, jitter check 6%, repeat search 5%. Nothing left is worth more than a
few percent without faster storage.

## 4. Timeline

| Commit | Content |
|---|---|
| `a0b24b3` | First complete tumor/normal Illumina pipeline |
| `340254e` | TINC, commons-math log-gamma, test scaffold |
| `0df08ef`, `d64de14` | Java unit tests ported |
| `08901b1` | High-depth and germline modes, panel-only fix, faster BAM decoding |
| `2f2d1dd` | SBX |
| `0b070c9`, `0eecd82` | Ultima, verified on public HG008 data |
| `3209491` | Append, split reads, MSI jitter parameters, multi-sample |
| `19b7fe2` | CRAM input; second optimisation round |

## 5. Java findings

Numbered S1 to S60 and grouped by kind. File references are to
`sage/src/main/java/com/hartwig/hmftools/sage/` unless they name `hmf-common`. "Reproduced" means the
port behaves the same way under `-java_compat`, so that its output equals Java's.

Since 2026-11-22 the port corrects most of these by default. A finding marked **Corrected by default**
takes the corrected path unless `-java_compat` is given; what the corrected behaviour is and what it
changes is in section 5.5. Every statement of sections 1, 2, 3 and 6 about equality with Java is about
the port run with `-java_compat`.

### 5.1 Wrong result

**Calling**

- **S1. `phase/PhasedGroupCollection.java:84`** `merge` uses `putAll`, which replaces this
  collection's list of groups for every key the other collection also has; those groups are lost.
  Reproduced. **Corrected by default.**
- **S2. `phase/PhasedGroupCollection.java:75`** `addPhaseVariants` tests whether the list for one key is
  empty, not the collection, so the first group under each later key replaces the collection's minimum
  position instead of lowering it. Collections that overlap can then stay apart. Reproduced.
  **Corrected by default.**
- **S3. `phase/PhasingUtils.java:175`, `:197`** `mergeUninformative` looks for the single variant of a
  group in `group.PositiveReadCounters`, the group itself, where the other groups were meant: the test
  is always true, and such a group is only removed when it is the only one. After
  `removeAll(matchingGroups)` the loop index is not adjusted, so the group that moves into the current
  place is skipped. Reproduced. **Corrected by default.**
- **S4. `filter/VariantFilters.java:474`** `belowMinTumorVafProbability` computes
  `supportCount * (int)round(min(altBQ / BQ, 1))`: the ratio becomes 0 or 1 before it multiplies. The
  comment above the line and the README give `SupportCount * min(AvgBQ[ALT] / AvgBQ[DP], 1)`.
  Reproduced. **Corrected by default.**
- **S5. `sync/FragmentSync.java:59-62`, `:106`** When `CombinedSyncData.formFragmentRead` throws, the
  cached first mate is not removed and the second mate is dropped. The first is processed alone at
  the end of the slice, out of order. Reproduced. **Corrected by default, intent uncertain.**
- **S6. `hmf-common` `region/BedFileReader.java:132-136`** `loadBedFileChrMap` starts a new list
  whenever the chromosome changes, so a BED file that comes back to a chromosome loses the earlier
  block. Reproduced. **Corrected by default.**
- **S7. `common/NumberEvents.java:36`** A CRAM normally leaves the `NM` tag out and htsjdk does not
  regenerate it, so `getOrCalcNm` counts the mismatches itself: only inside the region's reference
  bases, and for a merged fragment on the merged read where a BAM gives the first mate's tag. The
  count is not always the aligner's `NM`, so the event penalty on map quality changes for some reads
  and with it `RC_QUAL`, `AFMQ`, `MQF` and a few filters: CRAM and BAM of the same reads give
  different calls. On chr21 of the Illumina pair Java gives 6,023 records from the BAMs and 5,908 from
  the CRAMs. Reproduced. **Corrected by default.**
- **S8. `evidence/ReadContextCounter.java:629`, `:640`** With a split (RNA) read segment,
  `coversVariant` and `hasUncertainCoreBases` ignore the index they are given and use the segment's
  variant index, also when they are asked about the realigned index. Reproduced.
  **Corrected by default.**
- **S9. `evidence/Realignment.java:34`, `:57`** `checkRealignment` makes the realigned index relative
  to the split segment and then checks jitter with that index against the bases of the whole read.
  Reproduced. **Corrected by default.**

**TINC**

- **S10. `tinc/TincCalculator.java:37-54`** `calculate` has an early exit on `levelsSinceLowest`, which
  is never incremented; every test level is scored. Reproduced. **Corrected by default, intent
  uncertain.**
- **S11. `tinc/TincCalculator.java:130`** `populateDefaultLevels` accumulates the 0.05 steps in
  floating point without the rounding the two finer ranges get, so the levels are
  0.39999999999999997, ... and the last is 0.9500000000000003: a level of 1.0 is never tested.
  Reproduced. **Corrected by default.**
- **S12. `tinc/VariantData.java:116`** `isMnv` requires the alt to equal the ref, so it is never true:
  MNVs are looked up in gnomAD by their first base only and always pass the recheck of the germline
  alt support. Reproduced. **Corrected by default.**
- **S13. `tinc/VariantCache.java:128`** The TINC rewrite only queries chr1-22, X and Y, so the
  mitochondrial records are absent from the rewritten VCF. Reproduced. **Corrected by default.**
- **S14. `filter/FilterConfig.java:63`, `tinc/TincAnalyser.java:86-96`** With `-run_tinc` the first
  pass keeps filtered variants with up to 10 germline alt fragments rather than 3, and the rewrite
  drops those above 3 unless they are phased or were recovered. Reproduced in both modes: the two
  constants and `filterOutVariant` show this is the design (keep enough for the recovery, then return
  to the usual limit), not a fault. Left.

**Ultima**

- **S15. `quality/BqrCache.java:71`** `loadBqrFiles` overwrites the maximum raw quality for each sample
  while iterating a `HashMap`, so the value used is that of whichever sample id iterates last (and 0
  with `-skip_bqr`). Reproduced, bucket order included. **Corrected by default.**
- **S16. `evidence/ReadContextCounterFactory.java:77`** calls `markedInvalid()`, the getter, where
  `markInvalid()` was meant, so a context whose quality models fail to build is never marked invalid.
  Reproduced. **Corrected by default** (no effect on output: the counter counts nothing either way).
- **S17. `seqtech/UltimaRealignedQualModelBuilder.java:349`** passes the read index offset, not the
  read index, to `getStraddlingReadBases` for inserts. The bases are not read for an insert.
  Reproduced. **Corrected by default** (no effect).
- **S18. `hmf-common` `bam/ReadCigarState.java:34`** The copy constructor sets the valid flag to true
  whatever the source holds. Reproduced. **Corrected by default** (no effect: both callers copy a
  valid state).

**Append**

- **S19. `append/SageAppendApplication.java:135`** Existing records are fully decoded and written
  again by htsjdk, which re-formats them: `0.000` becomes `0.00`, `0.005` becomes `5.000e-03`,
  trailing missing genotype values are dropped, and standard header lines with a non-standard number
  or type are replaced. The appended sample keeps SAGE's own formatting, so one record holds both.
  Reproduced. **Corrected by default.**
- **S20. `append/SageAppendApplication.java:199`** Only the chromosomes that are processed are
  written, so the records of every other contig are dropped from the output. Reproduced.
  **Corrected by default.**
- **S21. `append/SageAppendApplication.java`** Append loads no Ensembl data, so the map quality cap of
  the highly polymorphic (HLA) genes, which the caller applies, never applies to the appended
  samples. Reproduced. **Corrected by default** when `-ensembl_data_dir` is given to `append`.

**Found by a review of the Java code (2026-11-25)**

S49 to S53 were found after the port was complete, by a review made for this purpose, and checked
against the code before they were taken in.

- **S49. `quality/MsiJitterCalcs.java:259-262`** `calcErrorRate` leaves out a change of more than 5
  repeat units with `impliedAltChange > MSI_JITTER_MAX_REPEAT_CHANGE`, but the change is negative for a
  deletion, so only insertions are left out. A deletion of 6 or more units of a microsatellite gets the
  model's floor (1e-4, Q40) and so the MSI quality path: no fixed indel penalty, the MSI quality used
  directly, no `minAvgBaseQual`, while the insertion of the same units gets the ordinary indel quality.
  The README ("indels of no more than 5 repeat units"), the model's range of -5 to +5
  (`MsiModelParams.java:56-59`) and `ReadContextCounter.java:204` (which takes `abs()` of the same value)
  all treat the two alike. Reproduced. **Corrected by default.**
- **S50. `dedup/DedupSnvMnv.java:22-54`** `dedupMnvOverlaps` takes the passing MNVs before its loop, but
  an MNV filtered by an earlier pair is still compared and can filter a third MNV that no passing MNV
  overlaps (README: "on PASS variants only"). Five or more substituted bases in one phase set are
  needed; the result is another representation (MNV against SNVs). Reproduced. **Corrected by default.**
- **S51. `tinc/VariantData.java:120`** TINC's `isLongIndel` compares the longer allele's length with 11,
  so an indel of 10 bases is long there, while the caller's test (`ReadContextCounter.java:255`) is the
  indel length. TINC then rechecks the germline filter of such an indel with other jitter counts than
  the caller set it with. Reproduced. **Corrected by default.**
- **S52. `sync/CombinedSyncData.java:126-129`, `common/NumberEvents.java:23-42`** A merged fragment
  copies every attribute of its first mate, NM included, and its event count is that NM less the
  indels of the merged cigar: the second mate's mismatches are left out, an indel of the second mate
  alone makes the count lower (even below zero), and the count depends on which mate starts first.
  Java's own path for a record without NM (a CRAM without the tag) counts the events on the merged
  read. Reproduced. Left: counting on the merged read is the obvious reading, but its effect is large
  (5.5) and the thresholds were set with Java's BAM behaviour; to be confirmed with the SAGE authors.

**Found in the upstream update (2026-11-27)**

- **S57. `evidence/Realignment.java:94` (from `de50386378`)** `realignedReadIndexPosition` moves the core end
  position back by the insert length for every insert whose read context's cigar ends in a soft clip. The move
  is right for a context built from a read whose insert is soft-clipped (no `I` for it, e.g. `14M84S` for a
  72-base insert on the public SBX slice): `ReadCigarInfo` counted the inserted bases as reference positions, so
  the core end position is the insert length too far. It is wrong when the context's cigar holds the insert as
  an `I` before the clip (e.g. `12M1I16M3S`): the position is then already right, the move puts the realigned
  index the insert length before the variant, and a read that needs realignment (aligned from the base after
  the insert, with the insert in a soft clip) loses its support. And it is too large by the length of the `I`
  when the context holds part of a long insert as an `I` and clips the rest (`27M1I21M51S` for 35 bases): there
  the extra base of shift makes the reads of a neighbouring allele one repeat base longer (a 36-base insert at the
  same `AAAT` repeat, called as its own record) match the context exactly and count as realigned support for
  both alleles. Unit tests: `tests::realign_soft_clip` (the three kinds of context). Effect: section 9. Reproduced.
  **Corrected by default**: the core end position is moved by the number of the variant's inserted bases that the
  context's cigar does not hold as an insert element (all of them, none, or the clipped part).

### 5.2 Nondeterminism

- **S22. `phase/PhaseSetCounter.java:12`** `LPS` ids come from a counter shared by all region threads,
  so a multi-threaded run does not reproduce itself (chr21, 8 threads against 1: 152 lines). The port
  assigns ids in region order in both modes, which is Java with `-threads 1`.
- **S23. `pipeline/RegionTask.java:235`, `:255`** `NEAR_HOTSPOT` is set after the region has been
  handed to the writer, so with several threads the flag can be missing. The port always sets it, in
  both modes.
- **S24. `pipeline/ChromosomePartition.java:42`, `:62`** The slice size is divided by
  `min(threads, 4)`, so a contig or a set of regions shorter than four slices (chrM) is cut
  differently at each thread count, and the calls change with it. Reproduced.
  **Corrected by default.**
- **S25. `sync/FragmentSync.java:43`, `phase/PhasingUtils.java` (`removeUninformativeLps`),
  `evidence/FragmentLengthCounts`** Leftover cached mates are processed, ties between phase sets are
  decided and fragment lengths are written in `HashMap` iteration order. The order is a function of
  the keys, so Java reproduces itself; the port emulates it in both modes. Left: no other order is
  more correct, and changing it would only move the result away from Java's.
- **S26. `SageCommon.java:33`** `setReadLength` iterates a `HashMap` keyed by chromosome, so the
  region it samples is arbitrary. The port samples the panel regions in chromosome order, in both
  modes.
- **S27. `dedup/IndelDeduper.java:293`** `checkDedupCombinations` enumerates a `HashSet` of objects
  without `hashCode`, so tied combinations are tried in identity-hash order and which of two variants
  is marked `dedupIndel` depends on JVM state. Observed: the same site gives opposite answers when
  the run's region starts 1.2 Mb earlier, or when a second reference sample is added. The port uses
  the group's order, in both modes.
- **S28. `filter/VariantFilters.java:415`** The binomial tail `1 - cdf` is exactly 0 in some cases only
  because of `FastMath` rounding, which selects a different penalty; one `MQF` value in a whole
  genome depends on it. The port evaluates the commons-math formulas with correctly rounded
  functions, which reproduces Java on every value checked. Left as it is in both modes.
- **S29. `vis/VariantVis.java:189`** The reads shown, and the order of reads that tie, come from
  `ThreadLocalRandom`: two runs never give the same file. The port shows reads spread evenly over the
  order they were read in, and ties keep that order, in both modes.
- **S30. `hmf-common` `vis/SvgUtil`** Text positions come from AWT font metrics, so the files depend on
  the fonts of the machine. The port has a table of the metrics Java reports on macOS, in both modes.
- **S31. `common/SageVariant.java:180`** The visualiser's "Filter" cell joins a `HashSet`; names that
  share a bucket come out in identity-hash order. The port writes those in the order the filters are
  declared, in both modes.

### 5.3 Crash

- **S32.** Java SAGE without `-sequencing_type` fails with a null pointer in the candidate stage on
  the public SBX chr21 slice (after "ignoring read with position before prior position"). The port
  runs through, in both modes.
- **S33.** Ultima reads without `tp` / `t0` tags: a null pointer in Java; the port stops with a message
  naming the read, in both modes.
- **S34.** The htsjdk in SAGE 5.1 rejects CRAM 3.1, which current samtools writes by default. The port
  reads 3.0 and 3.1, in both modes.
- **S35.** Visualiser: `-vis_variants` with a variant two or more regions out of position order stops
  the VCF writer ("Features added out of order"). The port does not stop, and its VCF is then out of
  order, in both modes.
- **S36. `vis/VariantVis.java:99`** Turning the visualiser on removes every chrM call: the chromosome
  length comes from `RefGenomeCoordinates`, which has none for chrM, the counters cannot be built and
  the VCF has no chrM records (0 against 128 on the Illumina pair), with nothing in the log. The port
  writes no visualiser files there and leaves the calls alone, in both modes.
- **S37.** With `-vis_purple_vcf` and no reportable variant, SAGE exits leaving an empty, unterminated
  `.vcf.gz`. The port writes the header and the index, in both modes.
- **S38. `vis/AminoAcidEvent.java:191`, `:253`** A stop written as `Ter` in a substitution
  (`p.Gly12Ter`) is looked up in the amino acid map and is a null pointer exception; an unknown
  three-letter code in a `delins` is appended as the four letters of `null`. Reproduced.
  **Corrected by default.**
- **S39. `vis/VisFileBuilder.java:243`** Any exception while a file is built means no file for that
  variant, with a warning. With `-vis_purple_vcf` that includes every reportable variant whose
  `IMPACT` transcript has no amino acid entry (intergenic and non-coding ones) or whose protein
  change the parser rejects. Reproduced. **Corrected by default.**
- **S54. `common/NumberEvents.java:36-44`** For a read without an `NM` tag, `getOrCalcNm` calls htsjdk's
  `SequenceUtil.calculateSamNmTag` on the region's reference bases (the region plus 1,000 bases), and
  htsjdk 2.23 indexes them without a bounds check. A split (RNA) read whose blocks reach beyond them,
  that is one with an intron longer than the margin, throws ("Exception counting mismatches for read
  ...") and the run stops: Java 5.1 cannot append (or call) a STAR BAM made without `NM`, which is
  STAR's default (`--outSAMattributes Standard`). Found on the tumor RNA of section 8 with its `NM`
  tags removed: Java exits on the first such read. The port used to count only the bases inside the
  region, in both modes; it now takes the bases of such a block from the genome, in both modes, which
  is what upstream SAGE does from `de50386378` on (5.1-beta.5). A read that Java can count is counted
  as before. On that input one record changes (`AFMQ` of the appended sample), and the port's
  differences between the BAM with and without `NM` are then exactly those of the beta.5 jar (4
  records, the merged fragments of S52). Upstream's fix covers the evidence stage only: the candidate
  stage (`candidate/RefContextConsumer.java:245`) passes a null genome, so the caller would stop with a
  `NullPointerException` on the same read; the port gives the candidate stage the genome as well. Unit
  tests: `tests::split_read::s54`. From the update to `883ab11bba` (section 9) the port counts as upstream's
  `NumberEvents.calcMutations` does (bases compared as bytes; a block whose bases cannot be had ends the
  count), so with `-java_compat` it matches Java there wherever Java does not stop.

### 5.4 Oddity

- **S53. `filter/VariantFilters.java:469-472`, `:293-299`** The "novel indel boost" of the VAF
  probability test adds 12 to the alt base quality, which only enters the alt to all quality ratio (of
  about 1 for indels), not the probability; and it applies to every Illumina tier. The README says that
  p is boosted, and only in hotspot and panel regions. Boosting p as the README says would pass novel
  indels at low VAF much more easily, so it is left as the code has it. Reproduced. Left.
- **S40. `vcf/VariantVCF.java:259`** `RC_JIT` is declared with three values in the header
  ("Shortened, Lengthened, QualityPenalty") and written with two. Reproduced.
  **Corrected by default.**
- **S41. `filter/SoftFilter.java:27`, `filter/VariantFilters.java:242`** The `fragmentLength` filter is
  declared in the header, but its test is commented out and no record can carry it. Reproduced.
  **Corrected by default** (the declaration is left out; no filter is invented).
- **S42. `evidence/JitterMatch.java:157-191`** `hasJitterMatchType` has a branch on
  `altMatchCount > 0`, which is never set. Dead code; left.
- **S43. `seqtech/IlluminaArtefactContext.java:35`** `requiresCheck` compares against `NO_BASE` (0)
  rather than `NO_INDEX` (-10). The value is not used; left.
- **S44. `hmf-common` `vis/BaseSeqViewModel.java:123-124`** `create` marks a base of a merged
  fragment as overlapped by looking both mates up at `baseIdx + unclippedStart`, where `baseIdx`
  counts read bases and not reference positions, so after an insert or a delete in the merged read
  the overlap border is drawn against the wrong columns. Reproduced. **Corrected by default.**
- **S45. `vis/VisFileBuilder.java:459-461`** The `AED` row is filled from the first counter only: its
  average edge distance over all reads under "Tumor" and over the alt reads under "Normal".
  Reproduced. **Corrected by default, intent uncertain.**
- **S46. `vis/VisFileBuilder.java:226`** `MQF` is shown as `round(round(x * 10) / 10)`, a whole
  number. Reproduced. **Corrected by default.**
- **S47. `vis/AminoAcidUtil.java:338`, `:405`** After a frameshift reaches a stop codon, an index
  into the fetched bases (`idx + 3`, `refIdx - 3`) is used as a genome position for the non-coding
  remainder of the track. Reproduced. **Corrected by default.**
- **S48. `vis/ReadTableColumn.java:43`** The reference names of a read and its mate are compared with
  `!=`. It works because htsjdk hands out one string per sequence; the port compares the indices.
  No effect.
- **S55. `evidence/ReadContextCounter.java:359` (`processRead`)** A read that reaches the variant only through
  a soft clip is matched against the read context; when its clipped bases support neither allele it is
  still counted in the depth (`DP`, `RC_CNT`), as a read that covers the variant. In RNA, STAR soft-clips
  read ends it cannot place, among them the overhangs of splice junctions, whose bases come from the next
  exon. On the 25 records of section 8 with RNA coverage, 178 such reads are counted in 18 records: 24
  support the alt, 1 the ref, 153 neither; leaving the 153 out would raise the RNA `AF` by more than 0.02
  on 5 records (at most 0.071). Only 60 of the 636 reads whose soft clip covers a variant end within 5
  bases of a splice junction seen in other reads, so junctions are not the main source. The code treats
  DNA the same way, and counting every read that covers the position may be the intent. Kept, in both
  modes; to be raised with the SAGE authors.
- **S56. `quality/QualityCalculator.java:136-150`** STAR gives a uniquely mapped read the map quality 255
  (the SAM value for "not available"). SAGE takes it as 255; the cap at 60 (`MAX_MAP_QUALITY`) applies
  only with the germline map quality ratio factor. The map quality side of a read's quality (255 less the
  event, improper pair and fixed penalties) then never binds for RNA reads, and `AFMQ` and `ARMQ` are near
  250. Taking 255 as 60 changes `RC_QUAL` on 4 of the 25 records of section 8 (by 0.2% to 0.7%) and
  rescales `AFMQ` on 18. Kept, in both modes: the effect is small, and appended samples are not filtered.

- **S58. `sage/README.md:70`, `:357`** The README lists `base_qual_fixed_penalty | 12` as an option, which the code
  has never read (it was a constant at `a224a2c` too), and gives `modifiedBaseQuality = baseQuality -
  baseQualityFixedPenalty (12) - edgeDistancePenalty`; from `de50386378` the penalty of an SNV or MNV is scaled
  from 6 to 12 by the recalibrated quality (12 from 36, 6 at 24 and below), and the fixed 12 holds only for
  indels, Ultima model qualities and Illumina artefact qualities. Documentation only; the port's README
  describes the scaled penalty. Left.
- **S59. `append/SageAppendApplication.java:280` (from `de50386378`)** `addAppendHeader` writes
  `##sageAppendVersion` only when the input has none, so a VCF appended twice keeps the first append's version.
  WISP (master) takes that one value to choose its quality-per-read threshold (21 from 5.1, 18 before) for every
  sample of the VCF, so samples appended by 5.1 after an earlier append by 5.0 are held to 18. With samples from
  two versions no single value is right (a version per sample would be needed). oncoanalyser appends once. Kept
  in both modes.
- **S60. `hmf-common` `variant/VcfFileReader.java:161` (from `de50386378`)** `parseToolVersion` reads a version
  as a decimal number ("5.1-beta.3" is 5.1, "5.10" is 5.1, so "5.9" is above "5.10"). Append's minimum (4.0)
  and WISP's 5.1 give the intended answer for every SAGE version so far and for a 5.10; a later minimum such as
  5.2 would not. Kept in both modes (`hmf_common::vcf::parse_tool_version`).

Also found while porting, in REDUX 2.1 (see `redux-rs`): `BaseQualityResults.addPerfCounter` merges a
counter into itself and fails during BQR; `-form_consensus` fails on the public Ultima slices on a
consensus read with an empty base.

### 5.5 Corrections

33 of the findings above are corrected by default. Each corrected place in the code asks
`compat::java("S<n>")` (`src/compat.rs`) and keeps the Java path beside the corrected one.

- `-java_compat`, or `SAGE_JAVA_COMPAT=1` in the environment: Java behaviour for every finding.
- `-java_compat_except S3,S7` (or `SAGE_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
  those listed. The effects below were measured this way, one finding at a time against
  `-java_compat`.

Both options are taken off the command line before the program parses it, and apply to the caller,
`sage-rs append` and `sage-rs tinc`.

The inputs are those of the regression suite: chr21 of the Illumina pair, and of the public SBX and
Ultima pairs, each with soft filters and with `-disable_soft_filter` ("unfiltered" below, where every
candidate is written: 110,762, 68,809 and 87,238 records against 6,023, 3,316 and 2,689). Effects are
counts of VCF records gone / new / changed. A record whose `LPS` ids are only renumbered, because a
phase set before it appeared or went, is not counted as changed; a record is counted when it shares a
phase set with other records than before, or when any other value differs (the read counts of its
phase sets, `LPS_RC`, included). "Intent uncertain" marks a correction where the Java
code allows more than one reading; the reading taken and the alternatives are given. These in
particular should be confirmed with the SAGE authors.

**Calling**

| Id | Corrected behaviour | Effect |
|---|---|---|
| S1 | on a merge the groups of both collections are kept; a group that matches one already there exactly is counted into it | with soft filters: Illumina and Ultima none; SBX 2 records new (kept because they are now phased with a passing variant), 9 changed, 1 phase set new. Unfiltered: 1,533, 3,682 and 1,536 records changed, of which 661, 1,832 and 666 share a set with other records than before; 12, 19 and 2 records move between PASS and a dedup filter (the dedup filters read the phase sets) |
| S2 | only the first group of a collection sets its minimum position | nearly the records of S1, since collections whose bounds are right merge before they can collide: with soft filters SBX 2 new and 9 changed; unfiltered 1,178, 3,428 and 1,667 changed (592, 1,735 and 673 in the membership of their sets); 12, 18 and 2 move between PASS and a dedup filter |
| S3 | a single-variant group is removed when no other group has its variant; the loop index follows the current group | the largest effect. With soft filters: 10 (Illumina), 0 (SBX) and 14 (Ultima) records changed, no FILTER changes. Unfiltered: phase sets 44,402 to 39,281, 31,301 to 29,927 and 35,387 to 33,463; 25,797, 18,577 and 15,962 records changed (20,237, 15,122 and 12,729 in the membership of their sets); 33, 2 and 4 records gone; 22, 9 and 12 move between PASS and a dedup filter |
| S4 | the product of support and quality ratio is rounded | `minTumorVAF` changes on 1, 3 and 1 records, none of them passing before or after; unfiltered none |
| S5 | **intent uncertain.** The pair is closed where the merge fails and the first mate is processed there, alone: Java's rule for every other pair that cannot be combined ("only the first record"), and the same reads as Java in effect counts. Alternatives: process both mates separately, as for pairs that do not overlap; or drop both | none: no merge fails on the six inputs |
| S6 | every block of a chromosome is kept, in position order | none possible with the resource files (no BED file comes back to a chromosome); unit test |
| S7 | htslib regenerates `NM` when it decodes a CRAM record | Illumina CRAM: 5,908 to 6,023 records (6 gone, 121 new, 2,905 changed in `RC_QUAL`, `AFMQ`, `MQF`; `minMQF` on 56, `maxGermlineRelQual` on 12, `minTumorQual` on 2; one more passing), now byte-identical to the BAM run; SBX CRAM none (already equal). CRAM 3.1 gives the same as 3.0 |
| S8 | with a split read, the index asked about is used, relative to the segment | RNA append (section 8): none on the slice of 25 regions; for the three whole tumor RNAs see section 8. The path is only reached by a split read that needs realignment |
| S9 | the jitter check of a split read is made on the segment | as S8; unit test |
| S24 | a contig is cut into the slices of one thread at any thread count | chrM at 3 threads: 2 of 128 records change in their read counts, and the output is byte for byte Java's at 1 thread; at 1 thread none |
| S49 | a change of more than 5 repeat units, longer or shorter, is outside the jitter model | Illumina chr21: 18 indel records changed (quality and support values), `FILTER` on 3, 245 passing as before. `-germline`: 1 gone, 335 changed (330 indels), `FILTER` on 22, passing 72,902 to 72,897. Whole genome, a 35x Illumina pair on Google Batch (2026-11-26; the image before and after S49 to S51, on the same BAMs): somatic 416,119 to 416,053 records (92 gone, 26 new, 7,274 changed), passing 4,378 to 4,365: 17 homopolymer deletions of 6 to 23 bases become `minTumorQual` with the same support, 2 of 6 bases become passing; after PAVE and PURPLE one of the 17 had been passing (3,824 to 3,823) |
| S50 | an MNV filtered on the way is skipped on both sides of a pair | Illumina chr21 `-germline`: `FILTER` on 4 records (2 MNVs, 2 SNVs), passing the same; with soft filters none |
| S57 | the core end position is moved by the variant's inserted bases that the read context's cigar does not hold as an insert element (upstream: by the insert length whenever the cigar ends in a soft clip) | somatic runs with soft filters (Illumina, SBX, Ultima chr21): none. Illumina chr21 `-germline`: 44 records (apart from `LPS` renumbering): realigned support up on 9 + 7 inserts, down on 2 + 3, one insert PASS to filtered; unfiltered SBX 14 records, PASS the same; unfiltered Ultima 2; append 1 to 17 records of the appended sample (section 9) |

**TINC**

Measured with the standalone analyser on the whole-genome VCF of the Illumina pair with simulated
contamination of 2%, 10% and 50% (`bench/simulate_tinc.py`; 352,620, 339,306 and 299,618 records in the
rewritten VCF), and on the two-tumor SBX runs.

| Id | Corrected behaviour | Effect |
|---|---|---|
| S10 | **intent uncertain.** The counter counts the levels since the lowest score and the search stops when five in a row have not improved on it. Alternatives: stop at the sixth (the test is `>=` on a counter the code never shows being incremented), or the full scan is what is wanted and the exit is dead code | none: the fitted level is the same at 2%, 10% and 50% (0.015, 0.100, 0.500) and on the SBX runs |
| S11 | the levels from 0.15 are 0.15, 0.20, ..., 1.00 | 2% and 10% none; 50%: the level is 0.5 in place of 0.49999999999999994, which moves the recovery limit: 1 record new, 2 change FILTER |
| S12 | an MNV is a variant whose ref and alt have the same length above 1: gnomAD is asked for each base and the lowest frequency taken (as `GnomadChrCache.getFrequency` does for an MNV), and the recheck of the germline alt support applies | simulations: 1, 12 and 71 records gone, 1,425, 957 and 275 changed (`GND_FREQ` on MNVs, FILTER), 2, 0 and 2 fewer recovered; SBX runs: 20 and 287 MNVs get `GND_FREQ` |
| S51 | an indel is long from 11 bases of indel length, as in the caller | Illumina chr21 with `-run_tinc` and the oncoanalyser options: none (no such indel is rechecked); unit test |
| S13 | records on contigs other than chr1-22, X and Y are written where the input has them; they take no part in the fit and have no gnomAD file; recovery and the germline alt support limit apply to them as they do in Java when the fit variants come from a file | simulations: 70, 75 and 93 mitochondrial records written, 2 of them passing; SBX runs none (chr21 only) |

**Ultima**

| Id | Corrected behaviour | Effect |
|---|---|---|
| S15 | the maximum raw quality is the maximum over the samples' recalibration tables (intent: the name; the alternative, one value per sample, would need the quality calculator changed). With `-skip_bqr` it stays 0: there is no table to take it from | none on the Ultima pair and on a run with two tumors: the tables of these samples have the same maximum |
| S16 | the context is marked invalid | none possible: the counter counts nothing either way |
| S17 | the read index is passed | none possible: the bases are not read for an insert |
| S18 | the copy keeps the valid flag | none possible: both callers copy a valid state |

**Append**

| Id | Corrected behaviour | Effect |
|---|---|---|
| S19 | the columns of the existing record are written as they were read; the appended samples follow in SAGE's formatting; a FORMAT field the record did not have (`LPSA`) is added at the end, with `.` or the default of `RegionAppendTask.checkGenotypeFields` for the existing samples; header lines are kept as read | formatting only: 1,793 of 3,316, 57,079 of 68,809 and 2,077 of 2,689 records differ, in `FSB` and `RSB` of the existing samples (`0.000` kept where Java writes `0.00`), `TQP`, and the place of `LPSA` in FORMAT |
| S20 | records of a contig that is not processed are written with an empty genotype for the appended samples, each contig where the input has it; records left out by `-specific_chr` stay out | none on the chr21 inputs; unit test |
| S21 | with `-ensembl_data_dir` the HLA genes are loaded and the map quality cap applies | none on chr21; appending a sample to the calls of the HLA region: 253 of 851 records change in `RC_QUAL` and `AFMQ` of the appended sample |

**Header and visualiser**

| Id | Corrected behaviour | Effect |
|---|---|---|
| S40 | `RC_JIT` is declared with two values | one header line |
| S41 | the `fragmentLength` filter is not declared | one header line less |
| S38 | `Ter` is the stop in a substitution; an unknown three-letter code is an error of the protein change | synthetic PURPLE VCF (119 variants): 4 variants with a stop written `Ter` get their file with the track; 5 with an unknown code no longer get a track with four extra amino acids (with S39 they get the file without a track) |
| S39 | when the amino acid track cannot be built the file is written without it, under the name a variant without a PURPLE record has, with a warning | synthetic PURPLE VCF: 119 files in place of 105; 173 flagged variants of the Illumina sample: 173 files in place of 155; the files Java does write are unchanged |
| S44 | the overlap of the mates is looked up at the column of the base | Illumina, 16 named variants: 8 files change; SBX and Ultima none (no mates) |
| S45 | **intent uncertain.** Each sample's column shows its own pair, "all reads, alt reads", as the `AED` tag of the VCF has it. Alternative: one value per sample (which of the two is not evident) | every file (16 of 16) |
| S46 | `MQF` is shown with one decimal | 11 of 16 files |
| S47 | the non-coding remainder after a frameshift's stop codon starts at the genome position after the codon | 4 of 105 files of the synthetic PURPLE VCF, 2 of 155 of the Illumina sample |

**All corrections together**, default against `-java_compat`:

| Input | Records | Gone / new / changed | Passing | Phase sets |
|---|---|---|---|---|
| Illumina chr21 | 6,023 to 6,023 | 0 / 0 / 11 | 245 to 245 | 77 to 78 |
| Illumina chr21, unfiltered | 110,762 to 110,730 | 33 / 1 / 26,564 | 103,358 to 103,344 | 44,402 to 39,369 |
| Illumina chr21, `-germline` | 110,833 to 110,833 | 0 / 0 / 16,677 | 72,902 to 72,900 | 28,888 to 26,250 |
| Illumina chr21, `-high_depth_mode` | 5,954 to 5,954 | 0 / 0 / 10 | 243 to 243 | 78 to 79 |
| Illumina chr21, tumor only | 110,833 to 110,833 | 0 / 0 / 16,677 | 72,933 to 72,931 | 28,905 to 26,265 |
| Illumina chrM (3 threads) | 128 to 128 | 0 / 0 / 10 | 2 to 2 | none |
| Illumina chr21 from CRAM | 5,908 to 6,023 | 6 / 121 / 2,905 | 244 to 245 | 76 to 78 |
| SBX chr21 | 3,316 to 3,318 | 0 / 2 / 12 | 174 to 174 | 55 to 56 |
| SBX chr21, unfiltered | 68,809 to 68,807 | 2 / 0 / 20,301 | 63,468 to 63,441 | 31,301 to 30,061 |
| Ultima chr21 | 2,689 to 2,689 | 0 / 0 / 15 | 150 to 150 | 49 to 51 |
| Ultima chr21, unfiltered | 87,238 to 87,234 | 4 / 0 / 16,816 | 83,344 to 83,343 | 35,387 to 33,470 |

With S49 to S51 (2026-11-25), Illumina chr21: 6,023 records, 19 changed (`FILTER` on 4), 245 passing;
`-germline`: 110,833 to 110,832 records, 358 changed (`FILTER` on 45), passing 72,902 to 72,895. The
output is the same at 1, 3 and 6 threads.

S52 is left as Java has it. Counting the events of a merged fragment on the merged read, its obvious
correction, changes Illumina chr21 from 6,023 records to 5,908 (121 gone, 6 new, 2,904 changed), 245
passing to 244, and with `-germline` from 110,833 records to 109,483 and 72,902 passing to 72,184: the
same as Java gives from a CRAM without NM tags (S7).
| SBX, two tumors | 92,402 to 92,402 | 0 / 0 / 15,111 | 68,212 to 68,211 | 22,826 to 22,179 |
| Ultima, two tumors | 10,715 to 10,716 | 0 / 1 / 17 | 150 to 150 | 42 to 43 |

The fragment length files (`-write_frag_lengths`) are unchanged on all six inputs that write them. In
the tumor/normal runs with soft filters no passing variant is lost or gained on any platform: the
changed records there are 1 to 3 with `minTumorVAF` added or removed on a variant that has other
filters (S4) and the rest in their phase sets. On chrM, `minTumorVAF` changes on 8 records, none
passing before or after. Where filtered candidates are written (unfiltered, `-germline`, tumor only)
the number of passing records moves by 2 to 27, through the dedup filters, which read the phase sets.

Whole genome (Illumina pair, `-run_tinc -write_tinc_vcf`, 6 threads, the same machine, one run after
the other):

| | Java mode | Corrected |
|---|---|---|
| Wall time, two runs each | 756 s, 744 s | 829 s, 766 s |
| CPU time (user) | 3,617 s, 3,600 s | 3,632 s, 3,594 s |
| Instructions | 41,163 G, 41,149 G | 41,136 G, 41,122 G |
| Peak memory | 0.84 GB, 0.71 GB | 0.77 GB, 0.76 GB |
| Records before TINC | 611,876 | 611,890 |
| Records in the rewritten VCF | 365,286 | 365,469 |
| Passing in the rewritten VCF | 15,882 | 15,884 |
| Phase sets | 4,419 | 4,352 |

The Java-mode output is byte-identical to the stored output of the port (both VCFs), and differs from
Java's own 16-thread run only in the 4,704 lines with thread-order `LPS` ids, as before.

Corrected against Java mode, before TINC: 14 records new, none gone, 1,000 changed (563 in the
membership of their phase sets; `minTumorVAF` added or removed on 269, none of them passing before
or after; the sample values of 4 records, all on chrM, through S24). The passing records are the
same 15,884. In the rewritten VCF: 191 new, of which 70 are the mitochondrial records (S13), 2 of
them passing variants that Java's rewrite loses, and 121 are filtered records now phased with a
passing variant, which the rewrite therefore keeps; 8 gone, filtered records that are no longer
phased and are above the rewrite's germline alt support limit (S14); 992 changed; `GND_FREQ` now on
318 MNVs (S12). The TINC level is 0 in both modes.

The corrected mode does the same work as Java mode: 0.07% fewer instructions, the same CPU time and
the same memory within the spread of two runs. Its wall times were the longer ones, but another job
was loading the machine and the disk (load average 10 to 15) and the two corrected runs differ from
each other by more than the modes do. The part of the run that is not parallel, the writer thread
that assigns phase sets, is not slower: with `-perf_summary` at 6 threads it is busy 0.8 s against
0.8 s on chr21, 6.4 s against 7.9 s on chr21 unfiltered and 3.7 s against 3.7 s on chr1, and the
elapsed times are equal (chr1: 56.7 s and 56.7 s corrected, 56.5 s and 64.3 s in Java mode). From a
CRAM the corrected mode needs about 9% fewer instructions, because htslib's `NM` replaces SAGE's own
count. All three corrected whole-genome runs gave the same two VCFs.

**CRAM against BAM.** Corrected, the Illumina, SBX and Ultima CRAMs give byte for byte the VCF of the
BAM of the same reads (CRAM 3.0 and 3.1 for the first two). In Java mode the Illumina CRAM differs
from its BAM in 121 + 6 records and 2,913 values, as Java does.

**Thread count.** In the corrected mode the output is the same at 1, 3 and 6 threads on nine runs
(Illumina chr21 with and without soft filters, chrM, SBX, Ultima, CRAM, TINC with two tumors, append).
In Java mode the same holds except for chrM (S24), where the port follows Java.

**Unit tests** of the corrected behaviour, each beside the Java behaviour: `phase::tests` (S1, S2,
S3), `tests::corrections` (S4, S5, S6, S15, S24, S38, S40, S41, S44, S49, S50), `tests::alignment_reading` (S7),
`tests::realignment` (S9), `tests::realign_soft_clip` (S57), `tinc::tests` (S10, S11, S12, S51), `append::tests` (S19, S20). Without one: S8
(the index of a split segment, covered through S9's test only for the jitter path), S13 and S21 (they
need input files), S16 to S18 (no observable effect), S39, S45, S46 and S47 (covered by the visualiser
runs above).

**Left as they are**, and why:

| Id | Why |
|---|---|
| S14 | not a fault: the higher limit with `-run_tinc` and the return to the usual one in the rewrite are both written out in the code |
| S25 | the `HashMap` orders are a function of the keys, so Java reproduces itself, and no other order is more correct; changing them would only move the result away from Java's. (S15, the one place where such an order picked a wrong value, is corrected) |
| S28 | the commons-math formulas are part of the result's definition; the port evaluates them with correctly rounded functions and agrees with Java on every value checked |
| S30 | the font table is what makes the files comparable with Java's at all |
| S42, S43 | dead code |
| S55, S56 | RNA in append: the intent may be what the code does (S55), the effect is small (S56); see 5.4 |
| S48 | no effect |
| S58 | documentation only |
| S59, S60 | upstream's behaviour; no effect on the cases of section 9 (one append per VCF; SAGE versions 4.2 and 5.1) |
| S22, S23, S26, S27, S29, S31 to S37, S54 | the port already does one fixed thing in both modes, because Java's own result is not reproducible or Java fails |

## 6. Later additions

- **`-write_frag_lengths`** (caller and append): `<output>.frag_lengths.tsv.gz` is byte-identical to
  Java's on SBX, Ultima, Illumina and two-tumor runs (1.1 to 5.3 million lines each). Java writes each
  variant's lengths in `HashMap` order, which the port reproduces.
- **`-specific_var_file`** (caller and append): identical VCFs. The caller matches on position and
  alleles, append on position only, as in Java.
- **TINC with two tumors**: the pre-TINC VCF, the rewritten VCF and the fit-variant file are identical.
  Java uses the first tumor and requires exactly one reference sample.

- **The visualiser** (`-vis_variants`, `-vis_pass_only`, `-vis_max_support_reads`, `-vis_output_dir`,
  `-vis_purple_vcf`): one gzipped HTML file per variant, with the read table, the info tables and, for
  the reportable variants of a PURPLE VCF, the amino acid track. `src/vis_svg.rs` writes SVG as
  JFreeSVG 5.0.5 does and holds the font metrics; `src/vis_html.rs` the CSS builder and the HTML that
  j2html renders; `src/vis_aa.rs` the amino acid track; `src/vis.rs` the rest. A counter collects reads
  only for a variant that is shown, so a run without the visualiser does the same work as before
  (SBX chr21, one thread: 37.7 s against 37.6 s, same VCF).

  Java's files are not reproducible (the reads shown and the order of tied reads are random) and
  depend on the machine's fonts, so the comparison sorts the read rows within each match type block
  and was made against Java on the same machine; everything else must be equal byte for byte. With
  `-vis_max_support_reads -1` both sides show every read:

  | Input | Run | Files | Result |
  |---|---|---|---|
  | SBX, public | 9 named variants: SNV, MNV, insert, delete; with and without soft filters | 9 | equal |
  | Ultima, public | 6 named variants | 6 | equal |
  | SBX, public | `-vis_pass_only` over 400 kb | 3 | equal |
  | SBX, public | two tumors and two references | 3 | equal |
  | SBX, public | tumor only, `-vis_output_dir` | 3 | equal |
  | SBX, public | CRAM input | 3 | equal |
  | SBX, public | synthetic PURPLE VCF: 119 coding variants given made-up protein changes of every kind | 105 | equal; the same 14 variants fail in both |
  | Illumina | 16 named variants (paired reads, merged fragments, mate columns, duplicate counts) | 16 | equal |
  | Illumina | `-vis_pass_only` over 2 Mb | 8 | equal |
  | Illumina | PURPLE VCF of the sample, and one with 173 of its variants flagged reportable | 4, 155 | equal; of the 173, the same 18 fail in both |

  With the default limits and with `-vis_max_support_reads 7` the two sides show different reads, so
  there the number of rows per block is compared, and everything outside the read rows: equal (SBX: 5
  files each; Illumina, default limits: 16; SBX `-vis_pass_only` over all of chr21: 174, which takes
  Java 2,392 s on one thread and the port 11 s on eight). The VCF of every one of these runs is
  identical to Java's.
  The unit tests take their expected values from Java (`src/tests/data/VisFixtures.java`): CSS, font
  metrics, SVG for base sequences, fragments, coordinates and gene tracks, 45 protein changes, and the
  amino acid track for 864 combinations of transcript, variant, protein change and view.

- **Base-quality bands** (2026-11-24; options of the port, not of Java): `-low_base_qual_threshold`
  and `-high_base_qual_threshold` replace the constants `LOW_BASE_QUAL_THRESHOLD` (26, candidate
  creation: `IlluminaUtils.lowQualInVariantCore`) and `ILLUMINA_HIGH_BASE_QUAL` (30, the strong
  high-quality support of `minTumorSupport` and the high-depth evidence test). Made for a whole-genome
  tumor whose recalibrated qualities stop below Q30, where every SNV candidate failed
  `minTumorSupport`. Checks, chr21 of the test pair: at the defaults,
  corrected mode is identical to the previous build; Java mode is identical to Java for the somatic,
  high-depth and germline runs; giving the defaults explicitly changes nothing; the same work
  (706.3 against 706.9 billion instructions). With the high band at Q60, which no base reaches, 10
  records pass instead of 245. Unit tests: both bands at candidate creation and in the support counts,
  and the option checks (0 to 93, low below high). The VCF header gains `##sageBaseQualThresholds`
  only when a band moves.

## 7. Not done

Untested, with what each test needs. None of these is known to fail; they have not been run.

| Gap | Why it is open | What a test needs |
|---|---|---|
| Ultima reads formed by REDUX consensus | REDUX fails on the public Ultima test data (section 5), so only non-consensus Ultima reads were compared | An Ultima BAM that REDUX has deduplicated with consensus, or a REDUX release with the failure fixed; then the usual Java-at-one-thread byte comparison |
| Reference genome version 37 | Every comparison used version 38 | A 37 reference, a 37 BAM pair and the 37 resource files (hotspots, panel, high-confidence, coverage, PON); run both tools with `-ref_genome_version 37` and byte-compare. Covers the chromosome naming without `chr`, the 37 PAR and MT handling |
| A Linux run | All builds and runs were on macOS arm64 | Build with `cargo build --release` on x86-64 Linux, run `cargo test --release` and the regression suite, and byte-compare one region against Java. Covers libdeflate and floating point on another target |
| Visualiser with several samples | The visualiser was compared with Java on one tumor and on tumor plus normal, not on a run with several tumor or reference samples | A run with two or more `-tumor` or `-reference` samples and `-vis_variants`; compare the HTML files with the read rows sorted, as in section 6 |
| A second deep exome pair | Prepared, not run | Run both tools on it at 16 threads and compare with `bench/compare_vcf.py` |
| Corrected mode on the whole runs of the second private project | S8 and S9 were measured on its three tumor RNAs (section 8); the DNA runs were not repeated in the corrected mode | Run the corrected mode on its DNA pairs and count the records that change |
| The corrections against the authors' intent | The corrected output cannot be compared with Java | Confirmation of S5, S10 and S45 (intent uncertain) and of S3 (the largest effect) by the SAGE authors |

Not done for speed: the visualiser builds its files on the writer thread, which only matters with `-vis_variants` on many variants. See the end of section 3 for the rest.

## 8. Append on RNA, and the 5.1 pre-releases (2026-11-27)

The oncoanalyser dev-3.1.0 line runs Hartwig's image `sage:5.1-beta.3`. The first combined DNA and RNA
run of the pipeline appended a STAR-aligned tumor RNA (duplicates marked by Picard MarkDuplicates) to
the PURPLE somatic VCF, once with that image and once with the port in Java mode: 25 of 4,384 records
differed, all in the RNA sample (`RC_QUAL` on 20, `AFMQ` on 5, `AFBQ` on 4), and the image wrote the
header line `##sageAppendVersion=5.1-beta.3`, which the port does not.

**Versions.** Upstream's tags: `sage-v5.1-beta.3` is `0e9f6ab950` (2026-11-13, on the branch of AUS-459),
`sage-v5.1-beta.4` is `a224a2cd8e` (master, 2026-11-17: the port's reference and the jar it was checked
against), `sage-v5.1-beta.5` is `de50386378` (master, 2026-11-19: the AUS-459 branch squashed in). beta.3
is not an ancestor of beta.4; it is an earlier build of the changes beta.5 brought to master. Built from
source, each logs "Sage version 5.1" (the pom's version); the image logs "5.1-beta.3".

**Reproduced** on a slice of the RNA (the reads within 2 kb of the 25 records) with the pipeline step's
command at 2 threads (records elsewhere lose their RNA reads in a slice, so only the 25 are compared with
the pipeline):

| Run | The 25 records | Whole VCF | Time, peak memory |
|---|---|---|---|
| 5.1 jar (beta.4) | the port's pipeline values | | 30 s, 2.5 GB |
| port, `-java_compat` | the port's pipeline values | byte-identical to the 5.1 jar, at 1 and 2 threads | 1.1 s, 75 MB |
| image 5.1-beta.3 (docker, amd64 emulated) | the image's pipeline values | 25 records and the header line differ from the 5.1 jar | 54 s |
| jar built from `de50386378` (beta.5) | 20 as the image, 5 as the port | header `##sageAppendVersion=5.1` | 28 s, 2.5 GB |

**Cause.** Two changes of the AUS-459 branch. Each was added on its own to an experimental build of the
port (not kept), and the pair reproduces the image:

1. The base quality penalty of SNVs and MNVs is scaled from 6 to 12 over recalibrated qualities 24 to 36,
   where it was a fixed 12 (`QualityCalculator.getBaseQualPenalty`; beta.3 and beta.5): 20 records
   (`RC_QUAL` on 19, `AFBQ` on 4). This RNA's qualities are binned at 2, 11, 25 and 37: a Q25 base keeps
   18.5 in place of 13, a Q11 base 5 in place of 0, and Q37 is unchanged.
2. beta.3 only: `NumberEvents.calcAdjustedNumMutations` returns the `NM` tag as it is, so an indel counts
   for its full length in the event penalty. beta.4 and beta.5 subtract length - 1 per indel (beta.5:
   "BWA counts indels for their full length, whereas Sage wants a mutation of 1 per indel"): 5 records
   (`AFMQ` on 5, `RC_QUAL` on 1). This is a fault of the pre-release that upstream has fixed.

With both changes the experimental build equals the image on all 4,384 records; with the first alone it
equals the beta.5 jar on all 4,384. So the port is right for the version it reproduces: the image is an
earlier pre-release, and 5 of the 25 differences are its own fault.

**Whether the port should follow.** Done on 2026-11-27 (section 9): the port's reference is now upstream master
`883ab11bba`, which has beta.5's SAGE code. The scaled penalty and the `NM` routine are taken as upstream has them,
append checks the input version and writes `##sageAppendVersion` (now that its qualities are made with the scaled
penalty, so that WISP applies the threshold meant for them), and the realignment change is taken with a correction
(S57). beta.3's `NM` is not: a fault upstream fixed. The pipeline's Java SAGE caller is still the beta.3 image, so
its calls differ from the port's by beta.3's `NM` fault (and by S57 and the other corrections in the default mode).

**RNA checks.** Spliced reads, STAR's map qualities, Picard's duplicates and soft clips, on the slice and
on the three whole tumor RNAs of the second private project (section 2):

- Spliced reads. A read whose splice gap covers the variant is not counted (Java and port). The
  corrections S8 and S9 change nothing on the slice; on the three whole RNAs they change no somatic
  record, and 0, 1 and 5 of the 1,910 to 1,957 germline records (`RC_JIT` on one; `AD`, `DP` and the
  quality fields on two; the append phase set `LPSA` on the others; one fragment length line).
- `NM`. These four RNA BAMs carry `NM` (they have the tags of STAR's `--outSAMattributes All`), and on
  the slice every primary, non-duplicate read's `NM` equals SAGE's own count on this reference. STAR's default leaves `NM` out,
  and Java 5.1 then stops: S54, now handled as beta.5 does. Without `NM` the port's output differs from
  the run with `NM` in 4 records, as the beta.5 jar's does, through the merged fragments of S52; with
  fragment sync off the two are identical.
- Map quality. 255 (unique) is used as it is: S56. Reads at 0 (more than 4 loci) are not evidence outside
  hotspots (`EVIDENCE_MIN_MAP_QUAL` 1); reads at 1 and 3 (2 to 4 loci) count in `AD` and `DP` with a
  quality of 3 at most.
- Duplicates. Reads flagged by MarkDuplicates are dropped, as REDUX's are (76% of the slice's reads).
  Picard forms no consensus, so `UMI_CNT` counts every read as having none.
- Soft clips: S55.

The existing append regression cases (9 on the public SBX and Ultima slices, 6 on the three tumor RNAs)
are byte-identical to the stored Java output, VCFs and fragment lengths, in Java mode at 2 threads after
the change for S54, and so are the caller on chr21 of the Illumina pair from BAM and from CRAM (where
Java mode counts `NM` itself). In the default mode the 15 append cases are identical to the build before
the change.

## 9. Following upstream to `883ab11bba` (2026-11-27)

The port's Java reference moved from `a224a2cd8e` (`sage-v5.1-beta.4`) to upstream master `883ab11bba`. In between,
SAGE changed in one commit, `de50386378` (`sage-v5.1-beta.5`), and hmf-common in the code SAGE uses there and in
`76f6774e65`; REDUX did not change. Each change was read in full and decided on its own (the decision table, for
`../UPSTREAM.md`, is in the merge notes):

| Change | Decision |
|---|---|
| Base-quality penalty of SNVs and MNVs scaled from 6 (recalibrated quality 24 and below) to 12 (36 and above) | accepted |
| SAGE's own `NM` count for reads without the tag (`NumberEvents.calcMutations`) | accepted; fixes S54 in the evidence stage |
| `realignedReadIndexPosition`: core end moved back by the insert length for inserts whose context cigar ends in a soft clip | accepted with a correction: S57 |
| Append: input SAGE 4.0 or later (read as a number, S60); writes `##sageAppendVersion` (S59) | accepted |
| hmf-common `SpecificRegions.includePosition` with no filters | accepted; no port has the old behaviour (every caller checks `hasFilters()` first) |
| Refactors, constants, documentation; other programs' parts | not applicable |

**Java references** were remade with the jar at `883ab11bba` for every regression case, with the same commands,
`-threads 1`, beside the old ones. With `-java_compat` the port is byte-identical to them on all 47 cases (VCFs,
fragment lengths, the TINC VCF and fit variants), at 2 threads; chrM differs at 2 threads by design (S24: 2 of 177
records) and is identical at 1.

**Effect in the default mode** (the port before the update against after, 2 threads; "changed" leaves out records
whose only change is a renumbered `LPS` id):

| Input | Records | Gone / new / changed | Passing |
|---|---|---|---|
| Illumina chr21 | 6,023 to 6,385 | 2 / 364 / 3,513 | 245 to 246 |
| Illumina chr21, `-germline` | 110,832 to 111,612 | 1 / 781 / 82,772 | 72,895 to 72,914 |
| Illumina chr21, unfiltered | 110,729 to 111,504 | 1 / 776 / 82,799 | 103,341 to 104,038 |
| Illumina chr21, `-high_depth_mode` | 5,955 to 5,967 | 1 / 13 / 489 | 243 to 243 |
| SBX chr21 | 3,318 to 3,374 | 4 / 60 / 1,717 | 174 to 174 |
| SBX chr21, unfiltered | 68,807 to 68,889 | 1 / 83 / 50,925 | 63,440 to 63,500 |
| Ultima chr21 (and unfiltered) | 2,695 (87,240) | none | 150 (83,349) |
| Append, 9 public cases | as before | 0 / 0 / 0 to 50,633 (the appended sample's values) | as before |
| RNA append, 6 cases | as before | 0 / 0 / 20 to 915 (RC_QUAL, AFBQ of the RNA sample) | as before |

The scaled penalty raises the final quality of bases below a recalibrated Q36 by up to 6: more low-quality
candidates pass the hard filters and are written as filtered records, and RC_QUAL / AFBQ rise. The Illumina tumor's
main BQR bins recalibrate to 36 or more, so only its low bins change; SBX's simplex bins (about Q32) change most;
Ultima SNVs use the BQR entry of the maximum raw quality and do not change. Every appended VCF has
`##sageAppendVersion=5.1`. The port does the same work as before (instructions within 0.5%).

**Truth set.** GIAB's draft HG008-T somatic small variant benchmark (v0.4, GRCh38, `tumorvariants.vcf.gz` inside
`all.bed`: 179 variants on chr21). PASS calls before and after are the same on both public slices: SBX 141 true, 8
false, 38 missed; Ultima 103, 19, 76 (also with calls under AF 0.05 dropped). In the unfiltered SBX run, ranked by
QUAL after removing calls with normal support, the 45 new records are all false and rank above a few true ones (true
variants before the 1,000th false one: 127 to 122). So the scaled penalty is taken as upstream's model on its
rationale; on these data it changes no call the benchmark can judge.

**S57.** With upstream's rule against the corrected one, both with the scaled penalty: the somatic runs with soft
filters are the same (Illumina, SBX, Ultima chr21); Illumina chr21 `-germline` 44 records, realigned support up on
9 + 7 inserts (normal, tumor) and down on 2 + 3, one insert from PASS to filtered; unfiltered SBX 14 records; append
1 to 17. Where upstream's rule gives more support it is through reads of a neighbouring allele (a 36-base insert
counted for a 35-base one at the same repeat) or a move by one repeat unit of a duplication (5.1, S57).

**Time and memory** (same machine, shared; Java at 1 thread; the port in Java mode):

| Case | Java, 1 thread | port, 1 thread | port, 2 threads | Java / port (1 thread) |
|---|---|---|---|---|
| chr21 | 506.1 s, 3179 MB | 59.3 s, 197 MB | 30.9 s, 271 MB | 8.5x |
| chr21_germline | 549.9 s, 4511 MB | 58.9 s, 162 MB | 31.7 s, 275 MB | 9.3x |
| sbx_chr21 | 114.8 s, 2639 MB | 38.3 s, 108 MB | 22.6 s, 144 MB | 3.0x |
| ultima_chr21 | 275.0 s, 2577 MB | 111.3 s, 109 MB | 60.6 s, 177 MB | 2.5x |
| cram_chr21 | 768.9 s, 2933 MB | 109.7 s, 132 MB | 56.7 s, 163 MB | 7.0x |
| multi_check_ultima | 730.1 s, 2697 MB | 287.2 s, 187 MB | 147.2 s, 239 MB | 2.5x |
| tinc_multi | 261.4 s, 19488 MB | 63.6 s, 174 MB | 35.6 s, 206 MB | 4.1x |
| append_check_nosoft | 321.9 s, 6180 MB | 65.8 s, 891 MB | 35.0 s, 919 MB | 4.9x |
| rna_t1_somatic | 185.8 s, 4290 MB | 40.5 s, 207 MB | 24.4 s, 225 MB | 4.6x |

