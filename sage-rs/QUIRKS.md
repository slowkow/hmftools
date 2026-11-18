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

## Where the port differs

- `dedup/IndelDeduper.checkDedupCombinations` enumerates subsets of a `HashSet` whose order comes from
  identity hash codes, so ties between equal variants are resolved arbitrarily in Java. The port uses
  the group's order. No difference has been observed.
- The binomial CDF follows commons-math3 3.6.1 (`Beta.regularizedBeta`, continued fraction) except that
  `logBeta` is computed from `lgamma`. Results can differ in the last bits; no difference has reached
  the output so far.
- Reads with `N` in the CIGAR (RNA) are skipped: `evidence/SplitReadSegment` is not ported.
