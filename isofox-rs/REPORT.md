# ISOFOX in Rust: record of the port

This report covers how `isofox-rs` was checked against Java ISOFOX 2.1, the faults found in the Java code, the
versions of ISOFOX involved, and the timings. Java ISOFOX 2.1 is this repository's
`isofox/target/isofox-2.1-jar-with-dependencies.jar`, built from sources that follow upstream commit ca04e6c39d.
Finding ids: IF. Options, outputs and tests are in [README.md](README.md).

The port also reproduces ISOFOX 2.0.1, the version for STAR-aligned RNA (`-isofox_version 2.0.1`), against a jar built
from upstream tag `isofox-v2.0.1`. Section 10 covers that mode; its findings are IF14 to IF18 in section 3.

oncoanalyser runs two ISOFOX steps: ISOFOX_QUANTIFICATION (`-functions
'TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS;RETAINED_INTRONS'`) and NEO_ANNOTATE_FUSIONS (`-functions NEO_EPITOPES`).
Both are ported.

## 1. What is ported

- **Fragment lengths:** the sampled distribution, with the `FRAG_LENGTH` and `FRAG_LENGTH_BY_GENE` files.
- **Fragment allocation:** each gene collection's reads are paired into fragments and matched to transcripts, giving:
  - the category counts and unspliced genes;
  - novel splice junctions, with the alt-SJ cohort;
  - retained introns;
  - chimeric read groups.
- **Expression:** the expected-counts fit (expectation maximisation) with the fragment length adjustment, the GC
  adjustment and its re-fit, raw and adjusted TPMs, the cohort medians and percentiles, panel TPM normalisation, and
  the summary statistics.
- **Fusions:**
  - the fusion finder (split, supplementary, discordant and realignment-candidate fragments);
  - cross-gene-collection and inter-chromosomal read groups;
  - the hard filter;
  - the pass filter, with the known fusion data and a fusion cohort.
- **NEO_EPITOPES:** the RNA fragments supporting NEO's fusion and point-mutation neo-epitopes, and the depth at their
  ends.
- **Options:** `STATISTICS` and every option of the run, as in Java. The other write types are accepted and write
  nothing (section 8).

The command line is Java's, and the bioconda/Hartwig wrapper's JVM options are accepted (`isofox -Xmx... -sample ...`).
The port is about 11,300 lines in 20 files, plus 1,200 lines of tests; NOTES.md maps the files to the Java classes.

## 2. How it was checked

### Defined Java reference

The reference is the jar run with the local JDK as `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2 ... -threads 1`,
one Java run at a time. Two things make Java's output vary, and this reference fixes both:

- **Thread timing (IF7).** At 2 threads Java's chromosome threads write the gene collection, alt-SJ, retained intron
  and fusion files in the order they finish. Fusion ids are numbered as threads reach them, and which fusions meet in
  the inter-chromosomal stage depends on finish order too.
  - On the slice, Java at 2 threads reorders 5 files and renumbers every fusion.
  - Sorted, its files equal those at 1 thread, and the fusion files also equal them once the id columns are dropped.
- **Identity hashes.** Several HashMaps that reach the output are keyed by objects without a hashCode (RegionReadData,
  the FusionFragmentType enum). Their iteration order then depends on identity hashes. With `-XX:hashCode=2`
  every identity hash is 1: insertion order up to 10 entries, then a tree ordered by compareTo. The port reproduces
  that order (NOTES.md). Without the flag, see IF13.

The port runs with `-java_compat`. Every output file is compared byte for byte, along with the exit code. Time and
peak memory come from `/usr/bin/time -l` for both. The run, comparison and input-making scripts are kept outside this
repository, with the outputs.

### Inputs

- **The slice (private):** a slice of a real tumour RNA BAM: 2 x 151 bp, aligned by STAR 2.7.3a as oncoanalyser
  does, sorted, and duplicate-marked by GATK MarkDuplicates.
  - It holds all of chromosome 22 plus 50 kb around both breakends of 30 of the pipeline's passing fusions (10 DEL,
    8 DUP, 8 BND, 4 INV).
  - Size: 3.28 M records. Read groups are renamed, and mates outside the regions are absent.
- **Run with:** the pipeline's exact command line and resources: the expected counts and GC ratios, the gene
  distribution and alt-SJ cohort, the driver gene panel and known fusion data.
- **Synthetic BAMs (no patient data):** generated from the public GRCh38 reference on public Ensembl genes
  (BCR, ABL1, EWSR1, FLI1, MAPK1, PPM1F, TOP3B, CRKL, SNAP29), 5,557 records:
  - spliced fragments of canonical and other transcripts, intronic fragments, exon skipping and novel splice sites,
    retained introns;
  - duplicates, unmapped mates, secondary alignments, low map qualities with XA tags, consensus (CR) reads, adapter
    read-through, soft clips matching the next exon;
  - a BCR-ABL1 translocation with supplementary alignments and discordant pairs, a SNAP29-CRKL deletion with split
    reads, an EWSR1-FLI1 translocation, inversions;
  - reads in the IGL region and on an alt contig.

  An empty BAM with the same header covers empty input.
- **NEO_EPITOPES:** NEO_FINDER outputs of two DNA pairs (HG008, public; the private pair, renamed), plus a NEO-format
  file made from the slice's own RNA fusions. Neither DNA pair is checked to be the RNA's subject, so those cases only
  exercise the code.

### Results (`-java_compat`)

| Case | Threads | Java time, peak RSS | Rust time, peak RSS | Result |
|---|---|---|---|---|
| slice, pipeline command | 1 | 291.2 s, 6.08 GB | 15.6 s, 1.54 GB | 9 of 9 files identical |
| slice, pipeline command | 2 | 152.6 s, 6.19 GB | 10.9 s, 1.57 GB | 9 of 9 identical to Java at 1 thread (Java at 2 threads: IF7) |
| synthetic, pipeline command | 1 | 208.4 s, 5.28 GB | 9.3 s, 1.14 GB | 9 of 9 identical |
| synthetic, pipeline command | 2 | 115.1 s, 5.27 GB | 7.2 s, 1.17 GB | 9 of 9 identical to Java at 1 thread; Java at 2 threads reorders gene_collection (IF7) |
| synthetic, `-specific_regions` (3 regions) | 1 / 2 | 14.7 s, 4.09 GB | 4.9 s, 0.82 GB | 9 of 9 identical |
| synthetic, `-specific_chr` (2) | 1 / 2 | 22.8 s, 4.25 GB | 4.9 s, 0.83 GB | 9 of 9 identical |
| synthetic, TRANSCRIPT_COUNTS only | 1 / 2 | 204.7 s, 4.76 GB | 6.4 s, 1.02 GB | 4 of 4 identical |
| synthetic, FUSIONS only | 1 / 2 | 115.2 s, 2.18 GB | 1.8 s, 0.27 GB | 2 of 2 identical |
| synthetic, STATISTICS only | 1 / 2 | 92.0 s, 2.17 GB | 1.6 s, 0.25 GB | 1 of 1 identical |
| synthetic, ALT_SPLICE_JUNCTIONS and RETAINED_INTRONS | 1 / 2 | 116.3 s, 2.44 GB | 2.3 s, 0.37 GB | 3 of 3 identical |
| synthetic, no expected counts, GC ratios or cohort files | 1 / 2 | 192.6 s, 2.48 GB | 2.2 s, 0.34 GB | 9 of 9 identical |
| synthetic, options (see below) | 1 / 2 | 204.0 s, 4.85 GB | 6.7 s, 1.17 GB | 10 of 10 identical |
| synthetic, `-exp_rate_frag_lengths` weights, `-drop_dups`, `-gene_read_limit 300` | 1 / 2 | 208.1 s, 5.25 GB | 7.2 s, 1.16 GB | 9 of 9 identical |
| synthetic, fragment lengths not matching the expected counts file | 1 / 2 | 10.2 s, 4.03 GB | 4.2 s, 0.61 GB | exit 1 in both; 6 of 6 identical (empty) |
| empty BAM | 1 / 2 | 188.8 s, 5.51 GB | 7.9 s, 1.47 GB | 9 of 9 identical |
| NEO_EPITOPES, HG008 neo-epitopes | 1 | 1.6 s, 1.25 GB | 0.7 s, 0.14 GB | identical, 27 rows (no RNA support: another subject) |
| NEO_EPITOPES, private pair's neo-epitopes | 1 | 1.7 s, 1.24 GB | 0.7 s, 0.14 GB | identical, 28 rows |
| NEO_EPITOPES, rows made at the slice's RNA fusions | 1 | 8.1 s, 2.15 GB | 2.1 s, 0.17 GB | identical: 119 rows, 28 with novel fragments |

- **Options case:** `-canonical_only -long_frag_limit 400 -fusion_min_frags_filter 1 -single_map_qual 50 -read_length 0
  -cancer_type Lung -frag_length_min_count 500 -write_type FRAG_LENGTH`.
- **Final binary:** every Rust figure above is from it. The Java figures are from the runs made during the work
  (the jar did not change).
- **Java's time on the small cases** is mostly loading: the expected counts file (107 MB) and the cohort files.
- **Load:** the machine (20 cores) was shared, with a load average of 7 to 8 from other work during the final runs.

### The default mode

Each case was also run with `-java_compat_except IFn` for each finding (on the slice: `finding_effects.sh`). The table
shows which findings change it. Every difference between the default mode and Java comes from these:

| Finding | Changes on the slice |
|---|---|
| IF2 | IntronicLength of every multi-exon gene: 25,629 gene_data rows |
| IF3 | 6 more alt-SJ rows, in 3 long genes (junctions more than 500 kb past the last gene start) |
| IF5 | 120 unfiltered alt-SJ rows become 123: junctions keep their own positions instead of being moved to an exon boundary |
| IF6 | 12 alt-SJ rows (the splice motif's strand now guides the choice of gene) |
| IF12 | the depth columns: 1,388 of 3,341 alt-SJ rows, 584 of 2,804 retained intron rows, 98 of 186 fusions. In pass fusions 23 lines change and one more fusion passes (56 against 55) |

IF1, IF4, IF8, IF9 and IF11 change neither the slice nor the synthetic cases. Unit tests show IF1, IF4 and IF9 in both
modes; IF8 and IF11 rest on reading the code. IF10 and IF13 keep Java's behaviour. The IF12 numbers are larger on the
slice than a whole BAM would give, because the slice's mates outside its regions are missing and more of its reads are
lone.

## 3. Findings

Each is corrected by default and reproduced by `-java_compat`, unless it says it is kept as Java.

- **IF1 (expected counts are scaled again when reused).**
  - What Java does: ExpectedCountsCache hands out its cached CategoryCountsData lists. ExpectedRatesData then scales
    them in place by the sample's fragment length frequencies (`CategoryCountsData.applyFrequencies`).
  - When a set is reused: a gene collection whose id does not match a set is given the first set containing all its
    genes, so a set can serve several collections. Typical causes are `-canonical_only`, `-specific_regions`, or
    Ensembl data newer than the file. A reused set is scaled again: frequencies squared, then cubed. With several
    threads, two can scale the same set at once.
  - Corrected: each use starts from the file's counts.
  - Effect: none in the cases (each set is used once). The unit test `expected_counts_frequencies_applied_once`
    shows the difference.
- **IF2 (IntronicLength off by the number of exonic regions).**
  - What Java does: gene_data's IntronicLength is the gene length minus the exonic bases. Both are measured as
    end - start, so the result is too large by the number of merged exonic regions minus one.
  - Corrected: inclusive lengths.
  - Effect: every multi-exon gene's IntronicLength. It is used by nothing else in the run.
- **IF3 (novel junctions far into long genes are dropped).**
  - What Java does: AltSpliceJunctionFinder accepts reads within 500 kb of the candidate genes, but takes the upper
    bound from the genes' latest start, not their end. Junctions more than 500 kb past the start of a gene longer than
    that are never evaluated.
  - Corrected: the genes' latest end.
  - Effect: 6 more rows on the slice.
- **IF4 (one missing transcript stops a chromosome's GC counts).**
  - What Java does: GcTranscriptCalculator.generateGcCountsFromFit warns about a transcript missing from the expected
    GC ratios file, then returns, which leaves the rest of the chromosome's gene collections out of the GC
    adjustment.
  - Corrected: only that transcript is skipped.
  - Effect: none with the shipped files (every transcript is present); it matters with Ensembl data newer than the
    file. Unit test `gc_counts_skip_only_the_missing_transcript`.
- **IF5 (homology test compares bases with themselves).**
  - What Java does: AltSpliceJunctionFinder moves a junction within 5 bases of an exon boundary onto the boundary
    when the bases moved match at both ends. Its test compares the start's bases with themselves, so every such
    junction is moved, homologous or not.
  - Corrected: the last (or first) d bases of the junction's two flanks are compared.
  - Effect: on the slice, 120 unfiltered alt-SJ rows become 123. A wrong move turns a novel junction into a known-boundary one,
    changing its type, its counts and its cohort frequency. Unit test `junction_homology_shift`.
- **IF6 (the splice motif never guides the gene choice).**
  - What Java does: AltSpliceJunction's donor-acceptor bases (GT-AG, CT-AC and the like) are never set. The motif's
    strand therefore never chooses between overlapping genes on opposite strands.
  - Corrected: set from the junction's base context.
  - Effect: 12 rows on the slice.
- **IF7 (output depends on thread timing).** Described under the defined reference above.
  - The port runs the chromosomes on `-threads` threads, but buffers their lines.
  - It runs the inter-chromosomal stage afterwards, in the single-thread order.
  - It numbers each chromosome's fusions after those of the chromosomes before it, and writes in that order.
  - So its output at any thread count equals Java's at one thread. Not switchable.
- **IF8 (the last fusion of a location is never compared).**
  - What Java does: FusionFinder.markRelatedFusions' inner loop runs `j < size - 1`, so the last candidate of a
    location is never compared. A location with exactly two candidates is not compared at all.
  - Corrected: `j < size`.
  - Effect: the RelatedSplicedIds/RelatedProxIds columns. None in the cases.
- **IF9 (discordant soft-clip junctions swapped).**
  - What Java does: FusionFragment.setNonJunctionData stores the soft-clip junction positions by junction side, then
    reads them by read index. With the reads in upper-lower order, the start and end junctions are swapped.
  - Corrected: read by side.
  - Effect: none in the cases. Unit test `dup_and_discordant_junction_order`.
- **IF10 (kept as Java: inversion homology).**
  - What Java does: FusionReadData.setJunctionHomology calls `reverseString(adjacentStr)` for same-orientation
    junctions and discards the result.
  - Why it is kept: a reversal would not be right either. The other breakend reads the opposite strand (a reverse
    complement), and the two breakends move in opposite directions, which the homology's use (merging close
    fusions) does not allow for. The intent is unclear, so Java's comparison is kept.
- **IF11 (one of the two discordant counts is lost).**
  - What Java does: FusionReadData.toFusionData sets DiscordantFrags from DISCORDANT and from DISCORDANT_JUNCTION
    fragments in turn, so the last one iterated wins. The order is identity-hashed (insertion order under
    hashCode=2).
  - Corrected: the sum, which also goes into the total fragments and the filters.
  - Effect: none in the cases. A fusion has both types only when some discordant pairs carry a soft clip at the
    junction.
- **IF12 (lone reads' depth counted twice).**
  - What Java does: FragmentAllocator adds the depth of reads whose mate was not found in the gene collection
    (processIncompleteReads). ChimericReadTracker.postProcessChimericReads adds the same reads again.
  - Corrected: once.
  - Effect: on a whole BAM, the reads of chimeric and cross-collection fragments. That raises the depth columns
    (alt-SJ, retained introns, fusions' coverage) and so the pass filter's allele fraction.
- **IF13 (kept as Java: output depends on identity hashes).**
  - Why Java varies: maps keyed by objects without a hashCode iterate in identity-hash order, and that order breaks
    ties. The keys are a read's exon regions (RegionReadData) and the fusion fragment types. The ties decide:
    - which of several equally good transcripts a junction or fusion reports;
    - which of two equally near exon boundaries a junction is moved to (IF5);
    - which region a retained intron's fragments are counted against;
    - which discordant count survives (IF11).
  - Evidence: the jar run twice on the slice without `-XX:hashCode=2` (the JVM's default identity hashes), at one
    thread. The two runs equal each other but not the hashCode=2 reference, in 4 files:
    - unfiltered alt-SJ: 529 rows differ only in the order of their transcript lists, and 36 in their values (the
      transcript or exon reported, and 2 junctions' positions);
    - retained introns: 2,805 rows instead of 2,804, with other fragment counts;
    - fusions: 16 rows differ in order only and 8 in values (the transcript or exon, and one fusion's discordant and
      total fragments);
    - pass fusions: 5 lines with another transcript or exon.

    The gene, transcript, gene collection and summary files did not change.
  - What else moves it: the identity-hash sequence depends on the JVM build, on how many objects were hashed
    before, and, with several threads, on which thread hashes them. So the pipeline's output (16 threads, the
    image's JVM) cannot be reproduced exactly by any run.
  - Why it is kept: the port gives the hashCode=2 result, which is one of the outputs Java can produce. A defined
    tie-break (for example by region start, then transcript id) would be a choice no Java run makes, and the intent
    is unclear. With IF11 corrected, the discordant count no longer depends on it.

The following findings are in ISOFOX 2.0.1 code that 2.1 does not share (section 10). They apply in 2.0.1 mode only.

- **IF14 (2.0.1: a split read whose junction is not found ends the run).**
  - What Java does: 2.0.1's ChimericUtils.findSplitReadJunction looks for the gap between consecutive mapped sections
    that equals the longest N. With a deletion beside the N (`10M5D100N10M`) the gap is wider and no junction is
    found. collectCandidateJunctions then dereferences the missing junction. The NullPointerException ends the
    chromosome task, and with it the run (exit 1).
  - Corrected: the read adds no candidate junction, as 2.1 does with its null check. (2.1 finds the junction from the
    cigar.)
  - Effect: only a split read that spans gene collections, or two genes' known splice sites, with a D beside its
    longest N. None on the slice or the synthetic cases. Unit test `split_junction_from_mapped_sections`.
- **IF15 (kept as Java: 2.0.1 fails on unpaired reads).**
  - What Java does: GeneRegionFilters.skipRead calls htsjdk's getMateUnmappedFlag for every read, which throws for
    an unpaired read. A single-end BAM fails at its first read.
  - Kept: 2.0.1 has no single-end path to correct to (2.1 added one). The port stops with an error naming the read.
- **IF16 (2.0.1: a fragment with two supplementary alignments loses two reads and never completes).**
  - What Java does: a supplementary alignment pairs with whichever read of the same name comes next (the fragment
    tracker is keyed by read name). When both reads are split, the fragment forms two pairs. addChimericReadPair
    puts the second in place of the first; the merge exists only under `-validate`, which logs it as an error.
    FusionReadGroup.isComplete is "three reads, or two without supplementary data", so four reads never complete.
  - Corrected: the second pair joins the group, and completeness counts the reads expected, as 2.1 does.
  - Effect: none on the slice. STAR's chimeric output there has at most one supplementary per fragment (16,627
    supplementary records, no read name twice). Unit test `fusion_group_completeness_201`.
- **IF17 (2.0.1: secondary alignments counted as unspliced or alternative fragments).**
  - What Java does: processIntronicReads counts a fully intronic pair as UNSPLICED (or ALT if split) even when both
    reads are secondary alignments. Every other path counts only non-secondary fragments, and TotalFragments excludes
    secondaries. Such a pair counts its fragment a second time.
  - Corrected: a secondary pair is not counted there. It still contributes to expression, with its map quality weight.
  - Effect on the slice: unspliced fragments 530,513 become 527,973, alternative 545,742 become 545,684 (204 gene
    collection rows), and UnsplicedFragmentPerc goes from 0.352 to 0.350.
- **IF18 (2.0.1: adapter trimming removes an aligned cigar element).**
  - What Java does: Read.trimAdapterSoftClipBases trims the smaller of the overhang past the mate and the soft clip.
    With no soft clip the trim is 0, but the last (or first) cigar element is still removed. Only the cigar string
    changes. Read.matches compares it when a read seen after one gene collection is met before the next, so such a
    read can be taken twice.
  - Corrected: nothing changes without a soft clip.
  - Effect: none on the slice or the synthetic cases. Unit test `adapter_trimming_201`.

IF1 to IF13 in 2.0.1 mode: IF1 to IF6, IF8, IF9 and IF11 are in 2.0.1's code too and are corrected the same way. IF10
and IF13 are kept, as in 2.1. IF7 is corrected the same way (2.0.1 has the same thread structure). IF12 does not arise:
2.0.1 has no processIncompleteReads, and counts a lone read's depth once.

### Kept as Java (oddities with no effect worth a finding, or choices)

- **RetainedIntronFinder** sets `spannedIsStart` before checking the candidate, so a rejected candidate can set it.
  The intent is unclear.
- **TranscriptExpression.calcTpmFactors** leaves out values below 0.5 only in the capped branch. The intent is
  unclear.
- **The expected counts header (unreachable):** with no configured fragment lengths, ExpectedCountsCache would take
  them from the file's header starting one column late. Configured lengths always exist, though: the option has a
  default, and an empty value fails to parse in both.
- **`-drop_dups`** is accepted and does nothing in 2.1, which always drops duplicates (section 4). In 2.0.1 it leaves
  duplicates out of expression, as there.
- **2.0.1's enriched genes** (RN7SL1-5 and RN7SK):
  - the genes, with 1 kb either side, become excluded regions, so their reads are skipped and the code meant to
    count them (processEnrichedRegionRead) is never reached;
  - EnrichedGenePercent is therefore the share of fragments in gene collections holding an enriched gene but outside
    that window;
  - those gene collections are left out of the GC ratio distribution, and the enriched genes out of the adjusted TPM
    total.

  The intent is unclear, so Java's behaviour is kept.
- **2.0.1 applies `-specific_regions` to mates and supplementaries** as well as to the read, so a pair with a mate
  outside the regions is skipped; 2.1 exempts them. The intent is unclear, so Java's behaviour is kept.
- **2.0.1 weights a low map quality fragment** in the category counts (0.5 to 0.1) but not in the GC ratio counts.
  The intent is unclear, so Java's behaviour is kept.
- **A failed run** (exit 1) leaves its output files empty, because Java exits without flushing its buffered writers.
  The port does the same: a writer not closed discards its buffer.

## 4. Versions: what the pipeline ran, 2.0.1 and 2.1

- **The pipeline did not run 2.0.1.** Its task log says `Isofox version 2.1-beta.5` and container
  `hartwigmedicalfoundation/isofox:2.1-beta.5`. The run's notes and the coordinator's status board say 2.0.1, which
  is wrong.
  - The ISOFOX sources at tag `isofox-v2.1-beta.5` equal this checkout's. The shared library differs only by a
    rename and an equivalent rewrite (EnsemblDataLoader, SpecificRegions).
  - The pipeline's image, run on the slice (1 thread, `-XX:hashCode=2`, linux/amd64 under emulation), writes 9 of
    9 files byte-identical to the jar's, and so to the port's.
- **ISOFOX 2.1 is not made for STAR BAMs.** Its README says: "To run Isofox on STAR-aligned BAMs, use Isofox v2.0 or
  earlier." 2.1 expects bwa-mem2 alignments against a transcriptome-augmented reference, lifted back by TARS and
  duplicate-marked by REDUX. What that means on the pipeline's STAR input:
  - **Duplicates are dropped** from every count: the BAM slicer discards them. Their number is meant to come back
    from REDUX's consensus reads (the CR tag), which GATK MarkDuplicates does not write. So TotalFragments excludes
    duplicates, and DuplicateFragments is 0 (the whole-sample summary shows 0). 2.0.1 counted duplicates for
    expression unless `-drop_dups` was given.
  - **Multi-mapped fragments are not recognised.** 2.1 finds them from bwa's XA tag, which STAR does not write.
    MultiMappedFragments is 0 throughout, and STAR's multi-mapping primaries (MAPQ 3 or less) count as unique, with
    full weight, in expression, novel junctions, retained introns and fusions. 2.0.1 took MAPQ 3 or less as
    multi-mapped and weighted it down.
  - **RN7SL1-5 and RN7SK are no longer excluded** (TARS unmaps them upstream). In this sample they come to about
    6,900 of the 1,000,000 raw TPM, so little changes.
- **2.0.1 on the slice** (`quay.io/biocontainers/hmftools-isofox:2.0.1--hdfd78af_0`, linux/amd64 under emulation,
  1 thread, the same command line):
  - **File names and headers:** 2.0.1 writes `retained_intron.csv`, not `.tsv` (the same columns), and names a
    column LowMapQual (gene collection) or LowMapQualFrags (gene and transcript data) where 2.1 says MultiMapped /
    MultiMappedFragments. Its summary also has EnrichedGenePercent. The other headers are the same.
  - **Counts (2.0.1 against 2.1):** TotalFragments 1,506,495 against 522,279. DuplicateFragments 984,216 against
    0; the difference is the duplicates, which 2.0.1 counts and 2.1 drops.
  - **Fragment shares:** spliced 0.256 against 0.330, unspliced 0.352 against 0.294, alt 0.362 against 0.338,
    chimeric 0.025 against 0.049. MedianGCRatio 0.475 against 0.485, spliced genes 606 against 593.
  - **Genes** (471 with 10 or more fragments in both): fragments per gene are a median 2.2 times higher in 2.0.1,
    and the Pearson correlation of log10(AdjTPM + 1) is 0.96. Of the 528 genes with a combined AdjTPM above 10, 115
    differ more than 2-fold. The slice's TPMs are normalised over its own genes, so the size of this differs on a
    whole sample.
  - **Other files:** unfiltered alt-SJ 1,847 rows against 3,341; retained introns 2,639 against 2,804; fusions 346
    against 186; passing fusions 46 against 55.
- **TARS and Ultima handling in 2.1:** 2.1 has no TARS or Ultima option or tag. Its support is the set of
  assumptions above: CR consensus counts, XA multi-mapping, and the enriched genes unmapped upstream. Ultima support
  is a single-end path for unpaired reads, which a paired STAR BAM never takes. The port has all of 2.1's paths.

## 5. Tests

`cargo test --release`: 65 tests, all passing. They are ports of the Java tests of the ported classes, on synthetic
genes and reads:

- ReadUtilsTest, ReadCountsTest, TransClassificationTest, FragmentTest;
- NovelJunctionsTest (alt-SJ types and retained introns), GcBiasTest, MiscTest, the ExpectationMaxFit case;
- FusionFragmentsTest (DEL and DUP fragments), FusionDataTest (fragment assignment to fusions).

Tests of the port's own:

- the findings in both modes: IF1 (`expected_counts_frequencies_applied_once`), IF2 (`intronic_length_finding`), IF4
  (`gc_counts_skip_only_the_missing_transcript`), IF5 (`junction_homology_shift`), IF9
  (`dup_and_discordant_junction_order`);
- the identity-hash region map (insertion order, then the compareTo tree);
- the SeqHashMap order against Java's, the parsers and the formatting;
- the 2.0.1 mode (`tests201.rs`, `version.rs`), each behaviour that differs from 2.1, in both modes where it applies:
  - map quality weights and multi-mapping, chimeric improper pairs, 2.0.1's adapter trimming (IF18), split
    junctions (IF14), realignment candidates;
  - the soft-clip rule (2.1's TransClassificationTest cases), excluded regions, the slice filters, write types;
  - fusion group completeness (IF16), the retained intron file and headers;
  - the allocator's counts and weights for duplicates, secondaries and low map quality;
  - the version option and the STAR header detection.

hmf-common-rs's suite (141 tests) passes with the additions, including a comparable-key tree order generated from a
JVM probe.

## 6. Speed and memory

Slice, pipeline command, final binary (Java: the jar, `-XX:hashCode=2`):

| Threads | Java | Rust | Rust is |
|---|---|---|---|
| 1 | 291.2 s, 6.08 GB | 15.6 s, 1.54 GB | 18.7 times faster, a quarter of the memory |
| 2 | 152.6 s, 6.19 GB | 10.9 s, 1.57 GB | 14.0 times faster, a quarter of the memory |

Rust's gain from the second thread is smaller than Java's: 4.7 s of its 15.6 s is loading the resource files, which
happens on one thread. On a whole sample the loading stays the same and the per-chromosome work, which runs in
parallel, grows about 45-fold.

The pipeline ran 2.1-beta.5 on the whole sample in 32 min with a 43.8 GB peak at 16 threads (Google Batch,
`-Xmx45g`). The whole-sample run of the port happens on Batch (not part of this work).

**Profile** (macOS `sample`, the slice at 1 thread). The run log's timestamps split the 15.4 s as follows:

| Phase | Time |
|---|---|
| Resource files | 4.7 s: the gene distribution 2.2 s (180 MB gzip), expected counts 1.2 s, alt-SJ cohort 0.8 s, GC ratios 0.5 s |
| Fragment length sampling | 1.1 s |
| Chromosome tasks | 9.1 s. Chromosome 22, the one fully covered, takes 3.0 s. Each of the other 23 takes 0.2 to 0.4 s, mostly building gene collections and reading their reference windows. The inter-chromosomal fusion stage takes under 0.1 s. |
| GC re-fit, TPMs and writing | 0.4 s |

- **By function** (shares of busy samples): BAM decompression (libdeflate) 14%, reference reads 7.6%, gzip
  decompression of the resource files 7%, the expectation-maximisation fit 5.3%. The rest is spread over text
  parsing, allocation, hashing and fragment allocation, with no single hotspot.
- **Earlier profiles** were led by the reference reads and by clearing the base-depth array. Three changes took
  the 1-thread run from 25.7 s to 16 s:
  - each chromosome's reference window is loaded when first used;
  - only the used part of the base-depth array is cleared;
  - the known fusion pairs are computed once per run.

**Memory.** The single-function synthetic runs show where the port's memory goes:
- the Ensembl data and the BAM: 0.25 GB (STATISTICS only);
- the alt-SJ cohort: 0.1 GB more;
- the expression resources (expected counts, GC ratios, gene distribution): about 0.8 GB more (TRANSCRIPT_COUNTS
  only, 1.02 GB).

The slice peaks at 1.54 GB. Java peaks at 2.2 to 6.2 GB with `-Xmx8g`; the pipeline gives it 45 GB on a whole
sample.

## 7. Shared code: additions to hmf-common-rs

- **`java_id_hash`:** `IdKey::compare_key` (default 0, so existing users are unchanged). A tree bucket of a key class
  that is Comparable in Java now orders by compareTo, with tieBreakOrder for ties, as `HashMap.treeify` and
  `putTreeVal` do. A test checks 40 key sequences against a JVM probe's output
  (`java_id_hash_comparable_expected.txt`).
- **New `java_seq_hash::SeqHashMap`:** a HashMap that iterates in Java's order for String and Integer keys, with
  constant-time insert, lookup and remove; the order is worked out when iterated. `clear()` keeps the capacity as
  Java's does. Tested against `JavaHashMap`. Used for the large, often-changed read-name maps.

No existing function changed behaviour.

## 8. Not ported

- **Write types:** EXON, SPLICE_JUNC, READ, FUSION_FRAGMENT, CHIMERIC_POSITION_DATA, SPLICE_SITE, TRANS_COMBO,
  GC_RATIO and MULTI_MAP_LOCI. These are debug and analysis outputs that oncoanalyser does not request; they are
  accepted and write nothing.
- **Cohort analyses** (`cohort.CohortAnalyser`, by `-analyses`):
  - SUMMARY;
  - splice junctions: ALT_SPLICE_JUNCTION, CANONICAL_SPLICE_JUNCTION, ALT_SPLICE_JUNCTION_MATRIX,
    SPLICE_VARIANT_MATCHING, RECURRENT_SPLICE_VARIANTS, SPLICE_SITE_PERCENTILES;
  - FUSION, PASSING_FUSION, RETAINED_INTRON;
  - expression: EXPRESSION_DISTRIBUTION, GENE_EXPRESSION_COMPARE, GENE_EXPRESSION_MATRIX, TRANSCRIPT_EXPRESSION_MATRIX,
    EXTERNAL_EXPRESSION_COMPARE;
  - PANEL_TPM_NORMALISATION (builds the panel normalisation file).
  - `fusion.cohort.FusionCompare` compares two fusion files.
- **Reference data:** `refdata.GenerateReferenceData` builds the expected counts (ExpectedCountsGenerator) and the
  expected GC ratios (ExpectedGcRatiosGenerator).
- **Panel normalisation** is ported: PanelTpmNormaliser, applied with `-panel_tpm_norm_file`.

## 9. Open points

- **The version and the aligner (section 4).** The pipeline ran ISOFOX 2.1 on STAR BAMs, which 2.1 says it does not
  support: duplicates and multi-mapped reads are treated as above. The port now reproduces both: 2.1 for bwa-mem3 +
  TARS/REDUX input, and 2.0.1 for STAR (section 10).
- **Adapter soft clips (2.0.1 and 2.1):** Read.likelyAdaperSoftClipping compares the signed insert size with the read
  length. A reverse mate's insert size is negative, so its soft clips always count as adapter, and never make a
  transcript ALT. Both versions do this; it is noted here rather than corrected, as a correction would change 2.1's
  default output as well.
- **The whole-sample comparison** with Java and the in-pipeline run happen on Google Batch.
  - Against the pipeline's own output (16 threads), expect two kinds of difference even with `-java_compat`: thread
    order (IF7: reordered rows and other fusion ids, equal once sorted and without the id columns) and identity
    hashes (IF13: tied transcript choices in alt-SJ, retained intron and fusion rows).
  - A byte-for-byte comparison needs Java run as the defined reference: `-threads 1 -XX:hashCode=2`.
- **Whole-sample time and memory** are not measured yet. The port's memory should stay near the slice's 1.5 GB plus
  the per-chromosome reads in flight.
- **Write types:** the debug and analysis write types are not ported (section 8). Port them if a downstream user needs
  one.
- **Inputs:** two items in the coordinator's ISOFOX baseline folder carry the RNA's real sample id in their names: a
  whole-BAM index beside the slice, and the raw pipeline outputs and logs. They were there before this work; rename
  or move them if every prepared file should say TUMOR.
- **IF10 and IF13** keep Java's behaviour because the intent is unclear. If Hartwig confirms the intent, the
  corrections are small: a reverse complement for inversion homology, and a defined tie-break for regions.

## 10. ISOFOX 2.0.1 mode

### Why

ISOFOX 2.1 expects RNA aligned by bwa-mem3, lifted over by TARS and duplicate-marked by REDUX. Hartwig's README asks
for ISOFOX 2.0 or earlier with STAR, and upstream oncoanalyser runs 2.0.1 after STAR. Section 4 has what 2.1 does to a
STAR BAM. The port therefore reproduces 2.0.1 as well. The reference is a jar built from upstream tag
`isofox-v2.0.1`; between it and 2.1 (ca04e6c39d) 103 files of ISOFOX and the shared library changed, mostly for TARS.

### Selection

`-isofox_version 2.0.1|2.1`, else `ISOFOX_VERSION`, else the BAM header (a STAR `@PG` line selects 2.0.1), else 2.1.
`-version` reports `2.0.1` or `2.1` accordingly (from the option or the variable). See README.md.

### What 2.0.1 does differently

All of the following is reproduced in 2.0.1 mode. The fragment allocator is a separate module (`allocator201.rs`);
the rest asks `version::v201()`.

- **Reads and fragments:**
  - **Duplicates** are read and counted in expression and the fragment counts, unless `-drop_dups` (or only fusions or
    novel junctions run). They add no depth, no novel junctions or retained introns, and no realignment candidates.
  - **Secondary alignments** are read when transcript counts run, each named `<read name>_<HI>`, so that a secondary
    pair forms its own fragment. They count in expression, weighted by map quality, but not in the fragment type
    counts, the strandedness or as chimeric fragments (IF17 is the exception).
  - **Multi-mapped** means a map quality of 3 or less (STAR's convention). Such a fragment has the expression weight
    0.5 (map quality 3), 0.33 (2), 0.2 (1) or 0.1 (0), and is counted as LowMapQual. The low map quality fragments
    are allocated to transcripts and genes as 2.1's multi-mapped ones are.
  - `-single_map_qual` defaults to 255, STAR's map quality for unique reads. It applies to the fragment length
    sampling and the neo-epitope reads, and to every read when transcript counts are not run.
  - **Supplementaries** pair with the next read of the same name (IF16).
  - **An improper pair** is chimeric.
  - **Regions and transcripts:** a read's exon regions and transcripts are worked out when it is read, before adapter
    trimming. The soft-clip rule checks the left side first, and has no exemption for short clips within the exon.
    Adapter trimming cuts a read's 3' soft clip back to its mate's unclipped 5' end (IF18).
- **Regions:**
  - `-excluded_regions`, and each enriched gene (RN7SL1-5, RN7SK by default) with 1 kb either side, are skipped:
    - the BAM is sliced around them, and a gene collection wholly within one is not read;
    - a read whose mate or supplementary lies in one is skipped.
  - `-specific_regions` filters mates and supplementaries too.
  - Gene collections with an enriched gene are left out of the GC ratio distribution, and the enriched genes out of the
    adjusted TPM total. The summary reports their collections' fragments as EnrichedGenePercent.
  - The fragment length sampling skips enriched genes, genes overlapping an excluded region, and genes with more than
    the maximum number of transcripts even with `-specific_regions`.
- **Chimeric reads and fusions:**
  - **Duplicates:** the known-splice-site test is skipped for duplicates. Duplicates join chimeric groups (without
    intronic transcript data or depth), and a fusion read group holding a duplicate is dropped once complete.
  - **Candidates and groups:**
    - a realignment candidate needs a soft clip of 3 to 10 bases;
    - a split read's junction comes from its mapped sections (IF14);
    - only a three-read group with a supplementary can be an inversion;
    - a fusion read group is complete with three reads, or two without supplementary data (IF16).
  - **The fragment builder:** it takes the first two reads with supplementary data (2.1 the read with the longer soft
    clip), never types a same-orientation pair as discordant within one gene collection, and makes realignment
    candidates of two-read fragments only.
  - **Junctions:**
    - the junction bases are 10 long either side (2.1: 20 before the junction), and the soft-clip window is at most
      10 bases;
    - there is no homology position adjustment and no adjustment towards known exon boundaries;
    - the junction types belong to the initial fragment: KNOWN from the exon validation, then CANONICAL from the
      splice motif, using both genes' strands by position when both genes are known;
    - the fusion keeps every transcript-exon of its initial fragment.
  - **Filters:** the hard filter counts every fragment type and has no supplementary map quality filter. The pass
    filter compares total fragments (2.1: split and discordant) with the requirement, and a short local fusion may have
    any orientation.
  - **Discordant fragments and fusion reads:** discordant fragments are always checked as DEL or DUP fragments, and a
    fusion read takes the exon references of its best match rank.
- **Novel junctions:**
  - a read is an alt-SJ candidate if it has more than one mapped section, inferred ones included;
  - its junction is the first gap between sections, whether a known splice site or not;
  - homology is tested only when the read meets no exon boundary;
  - retained introns are written comma-separated, as `retained_intron.csv`.
- **Outputs:** `LowMapQual` (gene collections) and `LowMapQualFrags` (gene and transcript data) for 2.1's multi-mapped
  columns, and the summary's EnrichedGenePercent before MedianGCRatio. The other files and columns are 2.1's.
- **NEO_EPITOPES:** 2.0.1 reads NEO's file by `-sample` and has no `-tumor` option. In 2.0.1 mode the port takes
  `-tumor` where given (as oncoanalyser passes it), else `-sample`.
- **Options:**
  - 2.0.1's: `-excluded_regions`, `-enriched_gene_ids`, and its write types;
  - not ported: READ_COUNTS, `-write_chimeric_reads`, `-write_chimeric_frags` (each ends the run with an error);
  - kept from 2.1: `-tumor` and `-log_read_ids`.

### How it was checked

- **Java runs:** the 2.0.1 jar with `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, one at a time.
- **The port:** `-java_compat -isofox_version 2.0.1`.
- **Comparison:** every output file byte for byte, and the exit code; time and peak RSS from `/usr/bin/time -l`.
- **Inputs:** the same private slice, synthetic BAMs and NEO_FINDER outputs as section 2, and the pipeline's command
  line (it passes no `-excluded_regions`). Two more cases:
  - the slice with upstream oncoanalyser 2.0.1's `-excluded_regions` (the resource file);
  - a synthetic case with excluded regions and enriched genes of its own.
- **NEO_EPITOPES:** Java 2.0.1 needs `-sample` to name NEO's file, so the cases pass the DNA sample id there. The
  port was also run with the pipeline's `-sample <rna id> -tumor <dna id>`, and its file compared with Java's.

### Results (`-java_compat`, 2.0.1 mode)

| Case | Threads (Java / Rust) | Java s, MB | Rust s, MB | Result |
|---|---|---|---|---|
| slice, pipeline command | 1 / 1 | 295.5, 6117 | 18.6, 1750 | 9 of 9 files identical |
| slice, pipeline command | 2 / 2 | 163.3, 6140 | 12.5, 1868 | 9 of 9 identical to Java at 1 thread; Java at 2 threads reorders 5 files (IF7: equal sorted, fusions also renumbered and 2 inter-chromosomal fusions with another initial read) |
| slice, plus upstream's `-excluded_regions` | 1 / 1 | 294.8, 6073 | 19.7, 1726 | 9 of 9 identical |
| slice, the jar against the first delegate's 2.0.1 image run | 1 | 295.5 | | 9 of 9 identical |
| synthetic, pipeline command | 1 / 1 | 202.9, 5171 | 8.8, 1131 | 9 of 9 identical |
| synthetic, pipeline command | 1 / 2 | | 7.0, 1162 | 9 of 9 identical to Java at 1 thread |
| synthetic, `-specific_regions` (3) | 1 / 2 | 14.0, 4482 | 4.7, 829 | 9 of 9 identical |
| synthetic, `-specific_chr` (2) | 1 / 2 | 21.9, 4470 | 4.7, 837 | 9 of 9 identical |
| synthetic, TRANSCRIPT_COUNTS only | 1 / 2 | 198.2, 4786 | 6.4, 1017 | 4 of 4 identical |
| synthetic, FUSIONS only | 1 / 2 | 100.4, 2240 | 1.8, 255 | 2 of 2 identical |
| synthetic, STATISTICS only | 1 / 2 | 88.7, 2163 | 1.6, 246 | 1 of 1 identical |
| synthetic, ALT_SPLICE_JUNCTIONS and RETAINED_INTRONS | 1 / 2 | 100.9, 2625 | 2.3, 381 | 3 of 3 identical |
| synthetic, no expected counts, GC ratios or cohorts | 1 / 2 | 189.9, 2499 | 2.2, 338 | 9 of 9 identical |
| synthetic, options (as the 2.1 case) | 1 / 2 | 201.7, 4924 | 6.6, 1198 | 10 of 10 identical |
| synthetic, fragment length weights, `-drop_dups`, `-gene_read_limit 300` | 1 / 2 | 201.6, 5171 | 6.9, 1163 | 9 of 9 identical |
| synthetic, fragment lengths not matching the expected counts | 1 / 2 | 9.8, 3507 | 4.1, 609 | exit 1 in both; 6 of 6 identical (empty) |
| synthetic, `-excluded_regions` (3 synthetic regions) and `-enriched_gene_ids` (2 genes) | 1 / 2 | 199.0, 5184 | 7.0, 1171 | 9 of 9 identical |
| empty BAM | 1 / 2 | 181.3, 5418 | 7.6, 1490 | 9 of 9 identical |
| synthetic BAM with a STAR `@PG` line, no `-isofox_version` | (Java: the synthetic case) / 1 | | | the port chose 2.0.1; 9 of 9 identical |
| NEO_EPITOPES, HG008 neo-epitopes | 1 / 1 | 1.5, 1232 | 0.7, 140 | identical, 27 rows; also with `-sample <rna id> -tumor <dna id>` |
| NEO_EPITOPES, private pair's neo-epitopes | 1 / 1 | 1.6, 1115 | 0.7, 140 | identical, 28 rows; also with `-tumor` |
| NEO_EPITOPES, rows made at the slice's RNA fusions | 1 / 1 | 8.0, 2256 | 2.2, 173 | identical, 119 rows (28 with novel fragments); also with `-tumor` |

- **Times:** the Rust times are from the final binary, on the shared machine (load average 4 to 6). Java's small
  cases are mostly the loading of the resource files.
- **Speed:** on the slice the port is 15.9 times faster than the jar at 1 thread and 13.0 times at 2, at 29% of its
  memory. It is slower than in 2.1 mode (15.6 s), because 2.0.1 also reads the duplicates (984,216 fragments on the
  slice) and the secondary alignments (268,158 records).
- **The pipeline's whole-sample 2.0.1 run** (Google Batch, 16 threads, the biocontainers image, 39 min 49 s, 43.2 GB
  peak) was checked against the port on the slice, where they share data: chromosome 22, which the slice holds
  whole.
  - The gene collection file's counts are equal for all 624 chromosome 22 collections, in every count column (total,
    duplicate, supporting, unspliced, alternative, chimeric, LowMapQual and both strands).
  - The 1,362 chromosome 22 unfiltered alt-SJ rows are the same junctions, with the same fragment counts and depths.
  - 11 alt-SJ rows report another tied transcript, and 125 and 140 rows list their transcripts in another order. The
    retained introns differ in 325 of 2,207 lines.
  - These are identity-hash ties (IF13). The pipeline's JVM runs without `-XX:hashCode=2`, and the jar run that way
    on the slice changes the same files by the same amounts (382 retained intron lines and 17 alt-SJ rows), and no
    others.

### The default mode in 2.0.1

Each correction alone, on the slice, against the jar:

| Finding | Changes on the slice |
|---|---|
| IF2 | IntronicLength of every multi-exon gene (gene_data) |
| IF3 | 6 more unfiltered alt-SJ rows |
| IF5 | 28 unfiltered alt-SJ lines change, and 2 more rows |
| IF6 | 12 alt-SJ rows |
| IF17 | 204 gene collection rows (unspliced 530,513 to 527,973, alternative 545,742 to 545,684), UnsplicedFragmentPerc 0.352 to 0.350 |

IF1, IF4, IF8, IF9, IF11, IF14, IF16 and IF18 change neither the slice nor the synthetic cases. IF14, IF16 and IF18
have unit tests in both modes.

### The 2.1 mode is unchanged

All of the first delegate's 2.1 cases were rerun with the final binary, in 2.1 mode (neither BAM names STAR): the
slice at 1 and 2 threads, the 13 synthetic cases and the 3 NEO_EPITOPES cases. All 18 are identical to the 2.1 jar's
outputs, as before.

### The version effect on the whole sample

The pipeline ran both versions on the same STAR BAM:

| | 2.0.1 | 2.1 (2.1-beta.5) |
|---|---|---|
| total fragments | 72.1 M | 22.9 M |
| duplicate fragments | 49.3 M (68%) | 0 |
| spliced / unspliced / alternative / chimeric shares | 0.234 / 0.370 / 0.371 / 0.019 | 0.317 / 0.312 / 0.340 / 0.045 |
| spliced genes | 20,778 | 20,607 |
| median GC ratio | 0.405 | 0.425 |
| unfiltered alt-SJ rows | 80,460 | 141,562 |
| retained intron rows | 124,664 | 127,083 |
| fusions | 38,074 | 20,837 |
| passing fusions | 629 | 391 |

### oncoanalyser

- **Selecting the mode:** for `--rna_aligner star`, the ISOFOX and NEO_ANNOTATE_FUSIONS steps should pass
  `-isofox_version 2.0.1`, and export `ISOFOX_VERSION=2.0.1` so that `isofox -version` (in versions.yml) reports it.
  A STAR BAM selects 2.0.1 by itself, but `-version` cannot see the BAM.
- **NEO_ANNOTATE_FUSIONS:** the step passes `-tumor <dna id>`, which Java 2.0.1 rejects (its parser stops at an
  unregistered option). Java 2.0.1 also reads NEO's file by `-sample`. The port accepts `-tumor` in 2.0.1 mode.
- **`-excluded_regions`:** upstream's 2.0.1 release passed `-excluded_regions`; the fork's ISOFOX module, written for
  2.1, does not. The port handles both.
