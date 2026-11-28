# Java behaviours found while porting, by subject

Each is a candidate for an upstream report. The ids (S1 to S60) are those of [REPORT.md](REPORT.md),
section 5, which has the Java file and line of each and, in 5.5, what the corrections change. File
references here are to `sage/src/main/java/com/hartwig/hmftools/sage/`.

Three states:

- **corrected**: the port takes a corrected path by default, and behaves as Java does with
  `-java_compat` (or `SAGE_JAVA_COMPAT=1`), so that its output equals Java's.
- **kept**: the port behaves as Java does in both modes.
- **port differs**: Java's result is not reproducible, or Java fails; the port does one fixed thing in
  both modes.

## Affect output

- **S57, corrected. Upstream's realigned index for inserts before a soft clip** (`evidence/Realignment`, from
  `de50386378`). For an insert whose read context's cigar ends in a soft clip, the core end position is moved
  back by the insert length. Right for a context built from a read whose insert is soft-clipped (the inserted
  bases were counted as reference positions); wrong for a context whose cigar has the insert as an `I` before
  the clip (the position was already right, and reads that need realignment lose their support); too large when
  part of a long insert is an `I` (reads of a neighbouring, one-base-longer allele in a repeat then count for
  both). Corrected: moved by the inserted bases the context's cigar does not hold as an insert.
- **S22, port differs. Phase set ids are thread-order dependent in Java.** `phase/PhaseSetCounter` is
  shared by all region threads, so `LPS` values (and anything ordered by them) change between runs with
  more than one thread. The port assigns ids in region order, which equals Java with `-threads 1`.
- **S23, port differs. `NEAR_HOTSPOT` is set after the region is handed to the writer**
  (`pipeline/RegionTask.finaliseResults`). With several threads the flag can be missing. The port always
  sets it.
- **S1, corrected. Lost phase groups.** `phase/PhasedGroupCollection.merge` uses `putAll`, which
  replaces this collection's groups for any key the other collection also has. Corrected: the groups of
  both are kept, and one that matches a group already there exactly is counted into it.
- **S2, corrected. Collection bounds.** `PhasedGroupCollection.addPhaseVariants` tests whether the list
  for one key is empty, not the collection, so the minimum position can be overwritten instead of
  reduced. Corrected: only the first group of a collection sets its minimum.
- **S3, corrected. Always-true test.** `phase/PhasingUtils.mergeUninformative` checks
  `group.PositiveReadCounters.contains(readCounter)` where a different group was presumably intended, so
  a single-variant group is only dropped when it is the only group. The loop index is also not adjusted
  after groups ahead of it are removed. Corrected: the other groups are searched, and the index follows
  the current group.
- **S4, corrected. Rounding before multiplying.** `filter/VariantFilters.belowMinTumorVafProbability`
  computes `supportCount * (int)round(min(altBQ / BQ, 1))`, so the ratio becomes 0 or 1 first.
  Corrected: the product is rounded.
- **S5, corrected, intent uncertain. Fragment sync exceptions.** When `sync/CombinedSyncData` throws,
  `FragmentSync` leaves the first mate cached and drops the second. The first is processed alone at the
  end of the slice. Corrected: the pair is closed where it fails and the first mate is processed there,
  alone, which is Java's rule for the other pairs it cannot combine.
- **S25, kept. Leftover cached mates are processed in `HashMap` iteration order**
  (`FragmentSync.emptyCachedReads`). The port reproduces Java's `String.hashCode` bucket order,
  including the table's peak capacity.
- **S25, kept. `HashMap<Integer>` iteration order decides ties** in
  `PhasingUtils.removeUninformativeLps`. Reproduced.
- **S6, corrected. BED loading starts a new list whenever the chromosome changes** (`hmf-common`
  `BedFileReader.loadBedFileChrMap`), so a BED file that revisits a chromosome loses the earlier block.
  Corrected: every block is kept, in position order.
- **S24, corrected. Partition size.** `pipeline/ChromosomePartition` divides the slice size by
  `min(threads, 4)`, so a contig shorter than four slices (chrM) is cut differently at each thread count
  and its calls change. Corrected: the slices of one thread at any thread count.
- **S26, port differs. Read length sampling** (`SageCommon.setReadLength`) iterates a `HashMap` keyed by
  chromosome, so the region it samples is arbitrary. The port samples panel regions in chromosome order.

## Do not affect output here, but look unintended

- S58, kept (documentation). SAGE's README lists a `base_qual_fixed_penalty` option (12) that the code never
  read, and its formula still has a fixed penalty of 12; from `de50386378` the SNV and MNV penalty is scaled from
  6 to 12 by the recalibrated quality.
- S59, kept. Append writes `##sageAppendVersion` only when the input has none, so after two appends the first
  version stays, and WISP chooses one quality-per-read threshold from it for every sample.
- S60, kept. hmf-common's `parseToolVersion` reads a version as a decimal number ("5.10" is 5.1).
- S40, corrected. `RC_JIT` is declared with three values in the VCF header and written with two.
  Corrected: declared with two.
- S41, corrected. The `fragmentLength` filter is declared in the header but its test is commented out.
  Corrected: not declared.
- S42, kept. `evidence/JitterMatch.hasJitterMatchType` has a branch on `altMatchCount > 0`, which is
  never set.
- S43, kept. `seqtech/IlluminaArtefactContext.requiresCheck` compares against `NO_BASE` (0) rather than
  `NO_INDEX` (-10); unused.

## TINC

- S10, corrected, intent uncertain. `TincCalculator.calculate` has an early exit that never triggers,
  because the counter it tests is never incremented; every test level is scored. Corrected: the search
  stops after five levels in a row that do not improve on the lowest score.
- S11, corrected. `populateDefaultLevels` accumulates 0.05 steps in floating point, so the last level is
  0.9500000000000003 and a level of 1.0 is never tested. Corrected: 0.15, 0.20, ..., 1.00.
- S12, corrected. `tinc.VariantData.isMnv` requires the alt to equal the ref, so it is never true: MNVs
  are looked up in gnomAD by their first position only and are never subject to the MNV germline alt
  support recheck. Corrected: the lengths are compared.
- S13, corrected. The TINC rewrite only queries chr1-22, X and Y, so mitochondrial variants are absent
  from the rewritten VCF. Corrected: they are written; they take no part in the fit.
- S14, kept (the design, not a fault). With `-run_tinc` the first pass keeps filtered variants with up
  to 10 germline alt fragments rather than 3, and the rewrite then drops those above 3 unless they are
  phased or were recovered.

## Ultima

- S15, corrected. `BqrCache.loadBqrFiles` overwrites the maximum raw quality for each sample while
  iterating a `HashMap`, so the value used is that of whichever sample id iterates last (and 0 with
  `-skip_bqr`). The port reproduces the bucket order under `-java_compat`. Corrected: the maximum over
  the samples (still 0 with `-skip_bqr`, where there is no table to take it from).
- S16, corrected (no effect on output). `ReadContextCounterFactory` calls `markedInvalid()` (the getter)
  where `markInvalid()` was meant, so a context whose quality models fail to build is never marked
  invalid; its counter has no matcher.
- S17, corrected (no effect). `createRealignedQualModels` passes the read index offset, not the read
  index, to `getStraddlingReadBases` for inserts. The bases are unused for inserts.
- S18, corrected (no effect). The `ReadCigarState` copy constructor resets the valid flag. Both callers
  copy a valid state.

## Split reads and append mode

- S8 and S9, corrected. With a split read segment, `coversVariant` ignores the index it is given and
  uses the segment's variant index, including for the realigned index; `checkRealignment` subtracts the
  segment start from the realigned index and then checks jitter with that index against the whole read.
  Quality scores and edge distance use the whole read, with whole-read indexes, in both modes.
  Corrected: the index given, relative to the segment; the jitter check on the segment.
- S19, corrected. Append rewrites existing records through htsjdk, which re-formats them: `0.000`
  becomes `0.00`, `0.005` becomes `5.000e-03`, trailing missing genotype values are dropped, and
  standard header lines with a non-standard number or type are replaced. The appended sample keeps
  SAGE's own formatting. Corrected: what the record held is written as it was read; a FORMAT field the
  record did not have (`LPSA`) is added at the end.
- S20 and S21, corrected. Append drops records on chromosomes it does not process, and loads no Ensembl
  data, so the HLA map quality cap never applies there. Corrected: those records are written with no
  values for the appended samples; with `-ensembl_data_dir` the HLA genes are loaded and the cap
  applies.
- S32, port differs. Java SAGE without `-sequencing_type` fails with a null pointer in the candidate
  stage on the full public SBX chr21 slice (after "ignoring read with position before prior position");
  the port runs through. Default-path comparisons on that data use chr21:20-40 Mb.
- **S54, port differs. Java 5.1 stops on an RNA BAM without `NM` tags.** For a read without `NM`,
  `common/NumberEvents.getOrCalcNm` counts mismatches with htsjdk's `calculateSamNmTag` against the
  region's reference bases (the region plus 1,000 bases), without a bounds check; a split read whose
  intron takes it beyond them throws and the run exits. STAR leaves `NM` out by default. The port takes
  the bases of such a block from the genome, in both modes, as upstream does from `de50386378` on; it
  used to count only the bases inside the region. Upstream's fix covers the evidence stage only (the
  candidate stage passes a null genome); the port gives both stages the genome. Since the update to
  `883ab11bba` the port counts as upstream's `calcMutations` does.
- **S55, kept, intent uncertain. Soft-clipped reads in the depth.** A read that reaches the variant only
  through a soft clip is counted in `DP` when its clipped bases support neither allele. In RNA these
  include the overhangs STAR clips at splice junctions; on the tumor RNA of REPORT.md section 8 they
  lower the RNA `AF` by up to 0.07. DNA is treated the same way.
- **S56, kept. STAR's map quality 255** (unique) is taken as 255, so the map quality side of a read's
  quality and its penalties never bind for RNA reads, and `AFMQ` / `ARMQ` are near 250. Small effect.

## CRAM

- S7, corrected. A CRAM normally leaves the `NM` tag out. htslib regenerates it when decoding; htsjdk
  does not, and `NumberEvents.getOrCalcNm` then counts the mismatches itself against the region's
  reference bases. Under `-java_compat` the port turns the regeneration off and counts as Java does; by
  default it lets htslib regenerate the tag, and CRAM and BAM of the same reads give the same calls. The
  count is not always the aligner's `NM` (which the BAM carries), so the event penalty on map quality
  changes for some reads and with it `RC_QUAL`, `AFMQ`, `MQF` and a few filters. On chr21 of the
  Illumina baseline pair, Java gives 6,023 records from the BAMs and 5,908 from CRAMs of the same reads;
  the port with `-java_compat` gives the same 6,023 and 5,908, and by default 6,023 from both. The SBX
  and Ultima slices give the same calls from BAM and CRAM.
- S34, port differs. The htsjdk in Java SAGE 5.1 rejects CRAM 3.1, the default of current samtools. The
  port reads 3.0 and 3.1 and gives the same calls from either.

## Visualiser

File references are to `vis/` and to `hmf-common` `common/vis/`.

Reproduced under `-java_compat`, corrected by default:

- S44. `BaseSeqViewModel.create` marks a base of a merged fragment as overlapped by looking both mates
  up at `baseIdx + unclippedStart`, where `baseIdx` counts read bases rather than reference positions,
  so after an insert or a delete in the merged read the overlap border is drawn against the wrong
  columns.
- S45, intent uncertain. `VisFileBuilder.renderSampleQualAndSiteInfo` fills the `AED` row from the first
  counter only: its average edge distance over all reads under "Tumor" and over alt reads under
  "Normal". Corrected: each sample's own pair, "all reads, alt reads".
- S46. `MQF` is shown as `round(round(mapQualFactor * 10) / 10)`, a whole number. Corrected: one
  decimal.
- S47. `AminoAcidUtil.getAltGeneRegionViewModels`, after a frameshift reaches a stop codon, uses an
  index into the fetched bases (`idx + 3`, `refIdx - 3`) as a genome position for the non-coding
  remainder.
- S38. `AminoAcidEvent.parse` appends the lookup result for each three-letter code of a `delins`, so an
  unknown code contributes the four letters of `null` as substituted amino acids. `p.Gly12Ter` (a stop
  written as `Ter` rather than `*`) and an insertion without a range are null pointer exceptions.
  Corrected: `Ter` is the stop; an unknown code and an insertion without a range are errors of the
  protein change, handled as in S39.
- S39. Any exception while a file is built means no file for that variant, with a warning. With
  `-vis_purple_vcf` that includes every reportable variant whose `IMPACT` transcript has no amino acid
  entry (intergenic and non-coding ones). Corrected: when the amino acid track cannot be built the file
  is written without it, under the name a variant without a PURPLE record has.
- S48, kept (no effect). The reference chromosome names of a read and its mate are compared with `!=`
  (`ReadTableColumn`, "TRL"). It works because htsjdk hands out one string per sequence; the port
  compares the indices.

Where Java's own output is not reproducible, and what the port does in both modes:

- **S29. Which reads are shown, and the order of ties.** `VariantVis.downsampleReadEvidenceRecords`
  draws the reads of each match type with `ThreadLocalRandom`, also when none has to be dropped, and the
  sort that follows is stable, so reads with the same orientation, combined quality and map quality come
  out in a different order on every run. The port shows all reads when they fit, otherwise ones spread
  evenly over the order they were read in (read `k * count / limit`), and ties keep that order. A
  comparison has to sort the read rows within each match type block, and can only compare the reads
  themselves when both sides show all of them (`-vis_max_support_reads -1` at a site with fewer than
  1000 reads per type).
- **S30. Text placement depends on the machine's fonts.** `SvgUtil` measures strings with AWT font
  metrics for the logical font "SansSerif". The port has a table of the advance widths, line height and
  descent that Java reports on macOS (where "SansSerif" is Lucida Grande), dumped from Java at the three
  sizes used. On a machine with other fonts Java's own files differ from these in the `x` and `y` of
  each `<text>`, the height of the coordinate header and the size of the box behind an indel length.
  Only printable ASCII is in the table; any other character is given the width of a space.
- **S31. Order of several filters.** The "Filter" cell joins a `HashSet<String>` collected from a
  `HashSet` of enum values. The names come out in hash-bucket order, which the port reproduces; names
  that share a bucket come out in identity-hash order of the enum values, which varies between JVM runs,
  and in the port in the order the filters are declared. The groups that share a bucket:
  `readStrandBias` and `realignFrequency`; `jitter` and `dedupMnv`; `maxGermlineAltSupport` and
  `minAvgT0Qual`; `minTumorQual`, `minAvgRawBaseQual` and `dedupMixedGermlineSomatic`; `minMQF` and
  `minGermlineDepth`; `fragmentStrandBias`, `dedupSnvMnv` and `dedupIndel`; `minAvgHpQual` and
  `minFragmentCoords`; `maxGermlineVAF` and `dedupMatch`.

Where the port differs, in both modes:

- **S36. Chromosomes outside 1-22, X and Y.** `VariantVis` takes the chromosome length from
  `RefGenomeCoordinates`, which has none for chrM, and the counter then cannot be built: with the
  visualiser on, Java writes no VCF records at all for chrM (0 against 128 on the Illumina pair) and
  logs nothing. The port writes no visualiser files there and leaves the calls alone.
- **S35. `-vis_variants` out of position order.** The regions are processed in the order the variants
  are given, so phase set ids follow that order. The VCF writer holds one region back and can place a
  region's records before those of the region just before it (the port does the same, and the VCFs are
  identical), but no further: when a variant is two or more regions out of place Java stops with
  "Features added out of order". The port does not stop, and its VCF is then out of order.
- `-vis_purple_vcf`: for a record with several alt alleles Java takes the one with the highest allele
  count; the port takes the first.
- S37. `-vis_purple_vcf` with no `REPORTED` record: both stop with "no reportable variants to
  visualise". Java leaves an empty, unterminated `.vcf.gz` without an index; the port writes the header
  and the index.
- Java loads every Ensembl transcript and amino acid sequence when the visualiser is on; the port loads
  those the reportable variants name. Files are built on the writer thread, after the region's records
  are final, where Java builds them on the region's thread.

## Where the port differs

In both modes.

- S27. `dedup/IndelDeduper.checkDedupCombinations` enumerates subsets of a `HashSet` of objects that
  have no `hashCode`, so the set's order comes from identity hash codes. Combinations that compare equal
  (same position, same ref length, same alt support) are then tried in that order, and the first that
  explains the indel decides which variants are marked `dedupIndel`. Java's answer therefore depends on
  JVM state: on chr21 of the public SBX data with two tumors and `-disable_soft_filter`, an SNV and an
  insertion at the same position swap `PASS` and `dedupIndel` when the same command is given a region
  that starts 1.2 Mb earlier, or a second reference sample. The port uses the group's order, which gives
  one of Java's two answers. Seen at one site in 92,054 records; a comparison that includes it can
  differ on those two lines.
- S28. The binomial distribution is a port of commons-math3 3.6.1 (`BinomialDistribution`,
  `SaddlePointExpansion`, `Beta`, `Gamma`), evaluated with correctly rounded `exp`, `log` and `log1p`
  (the `core-math` crate) in place of Java's `FastMath`. All 84 reference values generated by Java match
  bit for bit, but agreement is not guaranteed for every input. It matters because some results are
  differences from 1: `BinomialDistribution(54, 0.5).cumulativeProbability(0)` is exactly 0 in Java,
  which selects the fixed read strand bias penalty of 50 in `belowMinMapQualFactor`; the platform
  `exp`/`log` gave 1.1e-16 there and changed one `MQF` value in a whole genome.
- TINC only loads the gnomAD and PON entries at positions present in the VCF, in one pass per file; Java
  loads whole chromosomes. The output is the same.
- Ultima: the quality recalibration table is loaded once per process (Java merges into a static on each
  load); a failure to build quality models stops the run for a standard chromosome, where Java logs and
  continues with a counter that counts nothing; missing `tp` / `t0` tags stop the run with a message
  naming the read, where Java fails with a null pointer (S33).
- Append: header lines are kept as text; the GL to PL conversion, the INFO value-count check and the
  missing-header-line error of htsjdk are not reproduced.
- The tumor's reads are read from the file once per partition: the candidate stage keeps the records it
  reads, plus those within the evidence stage's 30-base soft-clip buffer, and the evidence stage takes
  its reads from there (Java queries the file again). A partition holding more than 32 MB of records
  (`-read_cache_mb`, per thread and tumor sample; 0 turns the cache off) is read a second time instead.
  The reads and their order are the same.
- Overlapping mates that are both a single match element are merged directly from the three runs (first
  read only, overlap, second read only) rather than by the base-by-base cigar walk; a unit test checks
  the two give the same read.
