# Java behaviours reproduced as-is

The port copies these on purpose so that output matches Java SAGE 5.1. Each is a candidate for an
upstream report rather than a silent fix. File references are to `sage/src/main/java/com/hartwig/hmftools/sage/`.

## Affect output

- **Phase set ids are thread-order dependent in Java.** `phase/PhaseSetCounter` is shared by all
  region threads, so `LPS` values (and anything ordered by them) change between runs with more than
  one thread. The port assigns ids in region order, which equals Java with `-threads 1`.
- **`NEAR_HOTSPOT` is set after the region is handed to the writer** (`pipeline/RegionTask.finaliseResults`).
  With several threads the flag can be missing. The port always sets it.
- **Lost phase groups.** `phase/PhasedGroupCollection.merge` uses `putAll`, which replaces this
  collection's groups for any key the other collection also has.
- **Collection bounds.** `PhasedGroupCollection.addPhaseVariants` tests whether the list for one key is
  empty, not the collection, so the minimum position can be overwritten instead of reduced.
- **Always-true test.** `phase/PhasingUtils.mergeUninformative` checks
  `group.PositiveReadCounters.contains(readCounter)` where a different group was presumably intended,
  so a single-variant group is only dropped when it is the only group. The loop index is also not
  adjusted after groups ahead of it are removed.
- **Rounding before multiplying.** `filter/VariantFilters.belowMinTumorVafProbability` computes
  `supportCount * (int)round(min(altBQ / BQ, 1))`, so the ratio becomes 0 or 1 first.
- **Fragment sync exceptions.** When `sync/CombinedSyncData` throws, `FragmentSync` leaves the first
  mate cached and drops the second. The first is processed alone at the end of the slice.
- **Leftover cached mates are processed in `HashMap` iteration order** (`FragmentSync.emptyCachedReads`).
  The port reproduces Java's `String.hashCode` bucket order, including the table's peak capacity.
- **`HashMap<Integer>` iteration order decides ties** in `PhasingUtils.removeUninformativeLps`. Reproduced.
- **BED loading starts a new list whenever the chromosome changes** (`hmf-common` `BedFileReader.loadBedFileChrMap`),
  so a BED file that revisits a chromosome loses the earlier block.
- **Read length sampling** (`SageCommon.setReadLength`) iterates a `HashMap` keyed by chromosome, so the
  region it samples is arbitrary. The port samples panel regions in chromosome order.

## Do not affect output here, but look unintended

- `RC_JIT` is declared with three values in the VCF header and written with two.
- The `fragmentLength` filter is declared in the header but its test is commented out.
- `evidence/JitterMatch.hasJitterMatchType` has a branch on `altMatchCount > 0`, which is never set.
- `seqtech/IlluminaArtefactContext.requiresCheck` compares against `NO_BASE` (0) rather than `NO_INDEX` (-10); unused.

## TINC

- `TincCalculator.calculate` has an early exit that never triggers, because the counter it tests is
  never incremented; every test level is scored.
- `populateDefaultLevels` accumulates 0.05 steps in floating point, so the last level is
  0.9500000000000003 and a level of 1.0 is never tested.
- `tinc.VariantData.isMnv` requires the alt to equal the ref, so it is never true: MNVs are looked up
  in gnomAD by their first position only and are never subject to the MNV germline alt support recheck.
- The TINC rewrite only queries chr1-22, X and Y, so mitochondrial variants are absent from the
  rewritten VCF.
- With `-run_tinc` the first pass keeps filtered variants with up to 10 germline alt fragments rather
  than 3, and the rewrite then drops those above 3 unless they are phased or were recovered.

## Ultima

- `BqrCache.loadBqrFiles` overwrites the maximum raw quality for each sample while iterating a
  `HashMap`, so the value used is that of whichever sample id iterates last (and 0 with `-skip_bqr`).
  The port reproduces the bucket order.
- `ReadContextCounterFactory` calls `markedInvalid()` (the getter) where `markInvalid()` was meant, so
  a context whose quality models fail to build is never marked invalid; its counter has no matcher.
- `createRealignedQualModels` passes the read index offset, not the read index, to
  `getStraddlingReadBases` for inserts. The bases are unused for inserts.
- The `ReadCigarState` copy constructor resets the valid flag.

## Split reads and append mode

- With a split read segment, `coversVariant` ignores the index it is given and uses the segment's
  variant index, including for the realigned index; `checkRealignment` subtracts the segment start from
  the realigned index and then checks jitter with that index against the whole read. Quality scores
  and edge distance use the whole read.
- Append rewrites existing records through htsjdk, which re-formats them: `0.000` becomes `0.00`,
  `0.005` becomes `5.000e-03`, trailing missing genotype values are dropped, and standard header lines
  with a non-standard number or type are replaced. The appended sample keeps SAGE's own formatting.
- Append drops records on chromosomes it does not process, and loads no Ensembl data, so the HLA map
  quality cap never applies there.
- Java SAGE without `-sequencing_type` fails with a null pointer in the candidate stage on the full
  public SBX chr21 slice (after "ignoring read with position before prior position"); the port runs
  through. Default-path comparisons on that data use chr21:20-40 Mb.

## CRAM

- A CRAM normally leaves the `NM` tag out. htslib regenerates it when decoding; htsjdk does not, and
  `NumberEvents.getOrCalcNm` then counts the mismatches itself against the region's reference bases.
  The port turns the regeneration off and counts as Java does. The count is not always the aligner's
  `NM` (which the BAM carries), so the event penalty on map quality changes for some reads and with it
  `RC_QUAL`, `AFMQ`, `MQF` and a few filters. On chr21 of the Illumina baseline pair, Java gives
  6,023 records from the BAMs and 5,908 from CRAMs of the same reads; the port gives the same 6,023
  and 5,908. The SBX and Ultima slices give the same calls from BAM and CRAM.
- The htsjdk in Java SAGE 5.1 rejects CRAM 3.1, the default of current samtools. The port reads 3.0
  and 3.1 and gives the same calls from either.

## Where the port differs

- `dedup/IndelDeduper.checkDedupCombinations` enumerates subsets of a `HashSet` whose order comes from
  identity hash codes, so ties between equal variants are resolved arbitrarily in Java. The port uses
  the group's order. No difference has been observed.
- The binomial distribution is a port of commons-math3 3.6.1 (`BinomialDistribution`,
  `SaddlePointExpansion`, `Beta`, `Gamma`), evaluated with correctly rounded `exp`, `log` and `log1p`
  (the `core-math` crate) in place of Java's `FastMath`. All 84 reference values generated by Java
  match bit for bit, but agreement is not guaranteed for every input. It matters because some results
  are differences from 1: `BinomialDistribution(54, 0.5).cumulativeProbability(0)` is exactly 0 in
  Java, which selects the fixed read strand bias penalty of 50 in `belowMinMapQualFactor`; the
  platform `exp`/`log` gave 1.1e-16 there and changed one `MQF` value in a whole genome.
- TINC only loads the gnomAD and PON entries at positions present in the VCF, in one pass per file;
  Java loads whole chromosomes. The output is the same.
- Ultima: the quality recalibration table is loaded once per process (Java merges into a static on
  each load); a failure to build quality models stops the run for a standard chromosome, where Java
  logs and continues with a counter that counts nothing; missing `tp` / `t0` tags stop the run with a
  message naming the read, where Java fails with a null pointer.
- Append: header lines are kept as text; the GL to PL conversion, the INFO value-count check and the
  missing-header-line error of htsjdk are not reproduced.
- The tumor's reads are read from the file once per partition: the candidate stage keeps the records
  it reads, plus those within the evidence stage's 30-base soft-clip buffer, and the evidence stage
  takes its reads from there (Java queries the file again). A partition holding more than 32 MB of
  records (`-read_cache_mb`, per thread and tumor sample; 0 turns the cache off) is read a second
  time instead. The reads and their order are the same.
- Overlapping mates that are both a single match element are merged directly from the three runs
  (first read only, overlap, second read only) rather than by the base-by-base cigar walk; a unit test
  checks the two give the same read.
