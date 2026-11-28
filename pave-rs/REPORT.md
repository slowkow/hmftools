# PAVE in Rust: record of the port

This report covers how `pave-rs` was checked against Java PAVE 2.0, the faults found in the Java code,
and the timings. Java PAVE 2.0 is this repository's `pave/target/pave-2.0-jar-with-dependencies.jar`,
built from sources that follow upstream commit ca04e6c39d of 2026-11-20; since 2026-11-28 the port follows upstream
master 883ab11bba (section 7), whose jar is the reference of `-java_compat`. Options and outputs are in
[README.md](README.md).

## 1. What is ported

All of `PaveApplication` as oncoanalyser runs it (PAVE_SOMATIC and PAVE_GERMLINE), with its options.

| Java | Rust |
|---|---|
| PaveApplication, PaveConfig, PaveConstants, FilterType | `main.rs`, `config.rs`, `constants.rs` |
| ChromosomeTask (task, applyFilters) | `task.rs`; reading the input in `input.rs` |
| VariantData, GeneDataCache, GeneCacheIndexing | `variant.rs`, `gene_cache.rs` |
| VcfWriter, TranscriptWriter | `output.rs` (the VCF writer itself is hmf-common-rs `vcf`) |
| impact: SpecificSpliceRegion and its resource `splice/special_splice_regions.tsv` (upstream 4adc014589) | `impact/special_splice.rs`, `resources/splice/` (built into the binary) |
| impact: ImpactClassifier, VariantTransImpact, CodingContext, CodingUtils, ProteinContext, ProteinUtils, HgvsCoding, HgvsProtein, SpliceClassifier, SpliceImpactType, PhasedVariants, PhasedVariantClassifier, VariantImpactBuilder, PaveUtils | `impact/` |
| annotation: ReferenceData, PonAnnotation, GnomadAnnotation, Mappability, MappabilityChrCache, ClinvarAnnotation, ClinvarChrCache, Blacklistings, Reportability | `annotation/` |
| hmf-common: PonCache, GnomadCache, VariantEffect, CodingEffect, VariantImpact(Serialiser), VariantTranscriptImpact, AltTranscriptReportableInfo, ReportablePredicate, Pathogenicity, PathogenicSummaryFactory, VariantTier, SpliceSites, TranscriptUtils.calcExonicCodingPhase, AminoAcids | in the modules above |

Not ported, as oncoanalyser does not run them:
- `pon_gen` (PonBuilder, the PON generator);
- `resources` (GnomadCacheBuilder, GnomadCacheCombiner, MappabilityBuilder).

PAVE draws no plots and calls no R or other programs.

### How the port reads

**Java** gives each chromosome task its own reader, which queries the input through its index. It loads
the resources a chromosome at a time and keeps all of that chromosome's PON, gnomAD and mappability
entries, so it peaks at 4-9 GB.

**The port** reads the input once:
1. It keeps the records the tasks would process (by filter type and specific regions).
2. It collects the positions they need, including each MNV base for gnomAD.
3. It reads each resource file once, in parallel with the others, keeping only the entries at those
   positions.

A lookup then gives what Java's gives:
- the PON, gnomAD and ClinVar lookups match exact positions (the first entry with the ref and alt, as
  Java's lists give);
- mappability takes the entry at or before the position, which a merge of the sorted positions with the
  file gives.

Java's ClinVar and mappability searches keep an index between variants. For positions in order (a
tabix-indexed VCF is sorted), they equal the stateless lookups.

### Java behaviour reproduced (not faults)

- **Output record order:** chromosomes in HumanChromosome order, and records in input order within each.
  The phasing never reorders them: a local phase set is flushed when the next variant's set differs.
- **Impact order:** impacts by gene in the Java HashMap order of the gene names (hmf-common-rs
  `JavaHashMap`), which orders PAVE_TI and decides ties in VariantImpactBuilder. Transcripts are in
  Ensembl cache order.
- **htsjdk 2.23's VCF writing:**
  - INFO keys and filters sorted;
  - genotypes re-encoded (VariantContextBuilder validates the genotypes, which decodes them);
  - QUAL as `%.2f` without `.00`;
  - doubles as formatVCFDouble;
  - header lines sorted, with the contig lines kept in their order;
  - a flag written as `KEY` (CLNSIG with an empty significance is written as a bare key, as Java does).
- **A comparison that never matches:** ImpactClassifier's `AltAminoAcids.equals(STOP_AMINO_ACID)`
  compares a String with a char. So a missense is always added for a protein change, and removed again
  once the stop gained is found. This has no effect, so it is not a finding.
- **"null" in strings:** RefGenomeSource returns null for an empty or inverted range, and Java's string
  concatenation appends "null". The port reproduces this where it can occur.

## 2. How it was checked

**Mode and threads:** every comparison ran the port with `-java_compat`. Both tools ran at 2 threads
unless stated, Java with `-Xmx8g`.

**What was compared:**
- the VCF after decompression;
- in every case below, the compressed BGZF bytes and the tabix index (decompressed and compressed) as
  well, and these were identical too;
- not `pave.version`, which holds the build time.

**The Java references** used the default JVM. Java's output does not depend on thread timing or identity
hashes: its maps are keyed by strings and used only for lookups, and the output is ordered by chromosome.

### Inputs

- **GIAB HG008** (public), from an oncoanalyser run of 2026-11-26 in which the Rust SAGE port wrote the
  VCFs:
  - somatic: 346,391 records, 17,254 passing;
  - germline.
- **A private whole-genome pair:** the SAGE VCFs of an all-Java oncoanalyser run, with the sample ids
  replaced by TUMOR and NORMAL in the file names and header.
- **Slices:** HG008 somatic chr17, chr21 and chr22 (bgzipped), and chr21 as a plain VCF.
- **Synthetic VCFs,** made from the public Ensembl cache, the reference and the resource files (no
  sample data). The main set has 365,463 records on chr5, chr17 and chr21.
  - **Where:** exon boundaries, splice donor and acceptor sites, coding starts and ends, UTRs, introns,
    upstream ranges and the TERT promoter.
  - **Kinds:** SNVs, MNVs, DELs and INSs up to 50 bases (with microhomology and repeats as SAGE annotates
    them), duplications, mixed INDELs, and local phase sets of nearby INDELs and SNVs.
  - **Resource entries:** the PON, gnomAD (and MNVs of neighbouring gnomAD SNVs), ClinVar (pathogenic
    and conflicting entries) and the blacklists.
  - **Tiers:** HOTSPOT, PANEL, HIGH and LOW, with some records filtered.

  A quarter of the main set (91,366 records) serves the option cases. Complex INDELs (both alleles
  several bases long, of other lengths) and longer ones are a separate set (PA4).

### Results

The transcript cases ran at 1 thread (see PA7). Time is wall time, and memory is the peak resident set
(`/usr/bin/time -l`). "Identical" means the decompressed VCF, the BGZF bytes and the tabix index
(decompressed and as written) all equal Java's, and so does the transcript TSV where one is written.

| Case | Java time, memory | Rust time, memory | Result |
|---|---|---|---|
| HG008 somatic, pipeline command | 21.4 s, 4.56 GB | 3.6 s, 0.20 GB | identical (17,254 records) |
| HG008 germline, pipeline command | 20.2 s, 4.12 GB | 2.5 s, 0.16 GB | identical (846) |
| private pair somatic, pipeline command | 20.4 s, 4.51 GB | 3.0 s, 0.17 GB | identical (4,206) |
| private pair germline, pipeline command | 19.1 s, 5.01 GB | 2.5 s, 0.16 GB | identical (805) |
| HG008 somatic, `-write_transcript_data` (Java 1 thread) | 25.2 s, 4.55 GB | 3.3 s, 0.25 GB | identical, TSV too |
| HG008 somatic chr17, chr21, chr22 | 18.5 s, 4.97 GB | 2.5 s, 0.16 GB | identical |
| HG008 somatic chr21 as plain VCF | exit 1 (PA3) | exit 1 | same failure; corrected: the chr21 records of the indexed run |
| HG008 germline, `-set_reportable` | exit 1 (PA9) | exit 1 | same failure and partial output; corrected: completes |
| synthetic, 365,463 records, somatic resources | 82.9 s, 8.98 GB | 9.5 s, 0.89 GB | identical (310,759) |
| quarter synthetic, `-filter_type ALL` | 80.4 s, 8.98 GB | 5.8 s, 0.37 GB | identical |
| `-filter_type PANEL` | 80.5 s, 9.01 GB | 5.5 s, 0.35 GB | identical |
| `-filter_type ALL -write_pass_only -set_reportable -write_detailed` | 80.9 s, 9.05 GB | 5.5 s, 0.37 GB | identical |
| `-write_transcript_data -set_reportable` (1 thread) | 109.3 s, 8.75 GB | 9.6 s, 1.14 GB | identical, TSV too |
| `-only_canonical` | 79.2 s, 9.07 GB | 4.4 s, 0.19 GB | identical |
| germline resources, `-filter_type ALL` | 19.7 s, 5.30 GB | 4.0 s, 0.36 GB | identical |
| `-specific_regions` (2 regions) | 23.9 s, 8.90 GB | 3.8 s, 0.17 GB | identical |
| `-specific_chr chr21;chr5` | 27.0 s, 9.02 GB | 4.5 s, 0.25 GB | identical |
| SBX PON, `-sequencing_type SBX` | 80.5 s, 9.06 GB | 5.6 s, 0.37 GB | identical |
| Ultima PON, `-sequencing_type ULTIMA` | 80.0 s, 9.05 GB | 5.6 s, 0.37 GB | identical |
| `-pon_artefact_file`, custom `-pon_filters` | 81.2 s, 8.98 GB | 5.3 s, 0.36 GB | identical |
| `-gnomad_freq_file` (one file, chromosome column) | 24.4 s, 6.81 GB | 3.2 s, 0.33 GB | identical |
| no annotation resources | 3.9 s, 1.69 GB | 2.0 s, 0.33 GB | identical |
| SBX-PON artefact INDELs (12,649 records, 1,806 with repeat count 7+) | 81.1 s, 8.97 GB | 3.9 s, 0.19 GB | identical |
| Ultima-PON artefact INDELs (16,778, 1,513) | 80.7 s, 8.99 GB | 3.9 s, 0.20 GB | identical |
| the SBX set run as ILLUMINA | 79.6 s, 8.94 GB | 3.9 s, 0.19 GB | identical |
| PA1 trigger (3 variants at a gene end), transcript TSV | 79.5 s, 8.91 GB | 3.7 s, 0.16 GB | identical, TSV too |
| complex INDELs (PA4, 9 records) | exit 1 | exit 1 | same failure; corrected: completes, 1 variant without impacts |
| Java PAVE's own output as input, `-set_reportable` | exit 1 (PA9) | exit 1 | same failure; corrected: completes |
| header only | 81.0 s, 8.96 GB | 3.6 s, 0.16 GB | identical |
| records on chrM and an unplaced contig as well | 80.4 s, 8.88 GB | 3.7 s, 0.16 GB | identical |

**Why the synthetic runs take Java 80 s, even with a header only:** this is Java's loading. Those inputs
carry no PON or gnomAD annotation, so Java reads the PON, the gnomAD file and the mappability of every
chromosome (3.5 s each), whether the chromosome has variants or not. SAGE 5.1's VCFs carry both (GND_FREQ
and PON_COUNT in their headers), so for the pipeline's inputs both tools take them from the records.

**The corrected (default) mode against Java,** records changed:

| Case | Records | Changed | By finding |
|---|---|---|---|
| HG008 somatic | 17,254 | 256 | PA2 256 (MAPPABILITY only), PA1 0 |
| HG008 germline | 846 | 2 | PA2 2 |
| private pair somatic | 4,206 | 90 | PA2 90 |
| private pair germline | 805 | 4 | PA2 4 |
| synthetic | 310,759 | 2,898 | PA1 1,091 (680 had no gene impact at all in Java), PA2 1,812, PA6 0 |

### Unit tests

`cargo test --release` runs 75 tests and ignores 25 (71 before upstream 883ab11bba; then SpliceContextTest's
testSpecialSpliceRegion, a test of the port that the special region stays with its transcript and exon, and two of
the resource reader). The Java unit tests are ported one module per class
(`src/tests/`): VariantTest, InsertCodingContextTest, DelCodingContextTest, MnvCodingContextTest,
NonCodingContextTest, SpliceContextTest, HgvsCodingTest, HgvsCodingAddressesTest, ProteinImpactTest,
PhasedVariantsTest, HomologyRealignTest and ResourcesTest.
- That is 85 tests, 25 of them ignored as in Java. One of HgvsCodingAddressesTest's ignored tests passes
  in the port and is left on.
- The helpers ImpactTestUtils, MockRefGenome and GeneTestUtils.createTransExons are ported too.

`src/tests/port.rs` tests:
- the corrected findings, in both modes;
- the PON and mappability readers;
- the gnomAD MNV frequency;
- the output record;
- the ClinVar pathogenicity.

In hmf-common-rs, `vcf` and `gz` have tests of their own: header order and repair, record parsing and
encoding, and the index read back through htslib.

### In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 8 on Google Batch ran PAVE_SOMATIC and PAVE_GERMLINE as the port, in Java mode, with every step
before them from the first run's cache. Against the first run's Java PAVE (Hartwig's 2.0-beta.2 image),
all 17,254 somatic and 846 germline records are identical. Only the `##PaveVersion` header line differs,
2.0 against 2.0-beta.2, and with it the compressed bytes and the index.

| Step | Port | Java |
|---|---|---|
| PAVE_SOMATIC | 6.0 s, 211 MB | 31 s, 6.3 GB |
| PAVE_GERMLINE | 4.8 s, 189 MB | 37 s, 9.8 GB |

## 3. Findings

Each is corrected by default behind `compat::java("PAn")`; `-java_compat` gives Java's behaviour.

- **PA1 (GeneDataCache.findGenes): a gene is dropped by the variant's end position.**
  - **What Java does:** with the task's indexing, a current gene is purged once its end (plus the
    upstream distance on the negative strand) is before the variant's END. The forward search then
    resumes past genes already seen.
  - **Effect:** a variant reaching past the end of a gene it overlaps (a DEL or MNV over a gene's last
    bases) loses the gene. So does every later variant in that gene (variants after a long DEL that
    starts in a gene and ends past it). The non-indexed search (no indexing) finds the gene.
  - **Evidence:** the PA1 trigger case, an SNV in a gene, then a DEL over its end, then an SNV before
    its end; Java writes no impact for the last two. Test `gene_search_keeps_genes_a_deletion_overlaps`.
  - **Count:** none on the four whole-genome inputs; 1,091 of 310,759 synthetic records gain impacts.
    Their INDELs reach up to 50 bases.
  - **Corrected:** genes are purged by the variant's start, before which no later variant lies.
- **PA2 (Mappability.loadEntries): the start is moved by 1 in a file that is already 1-based.**
  - **What Java does:** the mappability file lists the positions where the value changes, 1-based
    (`chr1 10001 10001 0.0023`; PAVE's MappabilityBuilder writes the 0-based start plus 1). The loader
    adds 1 again, so each entry becomes the inverted region (V+1, V), which contains no position.
  - **Effect:** a variant exactly at a listed position takes the previous entry's value. Between
    entries, the previous entry's value is meant, and taken.
  - **Count:** MAPPABILITY changes on 256 of 17,254 HG008 somatic records and 90 of 4,206 private
    somatic records, for example 0.061 against 0.047, or 0.333 against 0.077. PURPLE's somatic filtering
    reads it.
  - **Corrected:** positions as read. Test `mappability_positions`.
- **PA3 (ChromosomeTask with a plain VCF): every task reads the whole file.**
  - **What Java does:** a VCF whose name does not end in .vcf.gz is read without an index. Each
    chromosome task then iterates over all records, processing every human-chromosome record against its
    own chromosome's genes and resources.
  - **Effect:** the writer fails on the repeated records ("Features added out of order"), and Java exits
    1 with a truncated VCF.
  - **Evidence:** the chr21 plain VCF.
  - **Corrected:** each task takes its own chromosome's records, giving the indexed run's output. Java
    mode reproduces the failure (exit 1, the same error).
- **PA4 (classification): a variant the classifier cannot handle ends the run.**
  - **What Java does:** these throw:
    - complex INDELs, where both alleles are several bases long, of other lengths (ProteinUtils line
      113 indexes past the codon);
    - long DELs over several exons (CodingContext.trimToCoding gets an inverted range).

    ChromosomeTask turns the exception into `System.exit(1)`.
  - **Evidence:** found with gnomAD and ClinVar entries in the synthetic data. SAGE emits no such
    records in these inputs; it does emit simpler INDELs whose alleles differ in the first base, which
    work.
  - **Corrected:** the variant is written without transcript impacts, and a warning is logged.
- **PA5 (TranscriptWriter): an impact without a protein context gets one field for ten.** The
  placeholder `ProteinContext.empty()` is ",,,0,,,0-0,,,". It is written as one tab-separated field where
  the header has ten protein columns, so those lines are short of columns. Corrected: ten tab-separated
  fields with the same values.
- **PA6 (GnomadChrCache.getFrequency for MNVs): 0 means "no frequency yet".** The lowest of the bases'
  frequencies is kept with 0 as the start value, so a base of frequency 0 is replaced by the next one. No
  entry of the resource files is 0, so this has no effect with them. Corrected: the lowest.
- **PA7 (TranscriptWriter): the transcript file's line order depends on thread timing.** No switch.
  - **What Java does:** the chromosome tasks write their lines to the shared writer as they go. At 2
    threads, chr1's and chr2's lines interleave (603 runs on HG008 somatic) and the order differs between
    runs; the lines themselves are the same.
  - **The port** writes each chromosome's lines in chromosome order, which is Java's order at 1 thread.
    The Java references for the transcript file use 1 thread.
- **PA8 (ClinvarAnnotation.onChromosomeComplete): the chromosome caches are never freed.** No switch;
  memory only.
  - **What Java does:** the caches are keyed by the chromosome without "chr", but the removal asks for it
    with "chr". So all 2.2 million ClinVar entries stay in memory to the end.
  - **The port** keeps only the entries at the input's positions.
- **PA9 (annotations put without overwriting): an input that has the key already ends the run.**
  - **What Java does:** CLNSIG, CLNSIGCONF, BLACKLIST_BED, BLACKLIST_VCF, HOTSPOT and REPORTED are put
    with htsjdk's `putAttribute(key, value)`, which throws when the record already has the key.
  - **Effect:** two inputs end the run with exit 1:
    - an input annotated by PAVE before, at its ClinVar matches;
    - `-set_reportable` on any SAGE 5.1 record flagged HOTSPOT (2 HG008 somatic and 5 HG008 germline
      PASS records).
  - **Evidence:** HG008 germline with `-set_reportable`, and Java's own output as input.
  - **Corrected:** the value is replaced, as PAVE does for PAVE_TI and IMPACT.

Kept as Java does, as the intent is unclear:
- `-write_pass_only` ("Only annotate passing variants") drops the variants PAVE filters, not those the
  input filters: with `-filter_type ALL`, a record SAGE filtered is still written.
- When a variant hits two driver genes with other reportable transcripts, VariantImpactBuilder keeps
  those of the last such gene in HashMap order.
- **PA10 (upstream 4adc014589, SpliceClassifier.isWithinSpecialSpliceRegion): the special region's edges.** No
  switch; Java's behaviour in both modes, a question for the authors.
  - **Base -9 is in neither region.** PAVE's splice region covers intronic bases -1 to -8
    (SPLICE_REGION_INTRON_RANGE 8), the special MET exon 14 region -10 to -30, so an SNV at c.2888-9 is not in a
    splice region while those at -8 and -10 are. Upstream's test asserts it (-9 false), but as the region's edge,
    not as a choice about -9; the commit message says only "-30 -> -10". If the region is meant to extend the
    standard one, it would be -30 to -9.
  - **INDELs reach one base past the region.** The region is tested against the variant's Position to
    EndPosition, which for a DEL are the base before and the base after the deleted ones and for an INS the bases
    on either side, where the standard region tests the altered bases (altPositionsOverlap). So a DEL of -34 to -31
    (the base after it is -30) or of -9 alone (the base before it is -10) is in the region though no base of -30 to
    -10 changes; an SNV at -31 or -9 is not. Whether a deletion next to the region counts (it shifts the region
    towards the exon, for the -9 one) is a biological choice the code does not state.
  - **Evidence:** the synthetic MET case (`met14`, section 7): of 116 INDELs flagged, 4 alter no base of -30 to -10
    (c.2888-34_2888-31delCTTT, c.2888-32_2888-31delTT, c.2888-31delT, c.2888-9delT); the SNVs flagged are exactly
    -30 to -10.

## 4. Speed and memory

At 2 threads the port takes a sixth of Java's time on the pipeline inputs and about a twentieth of its
memory (table above).
- Java's time hardly depends on the number of variants: the synthetic header alone takes 81 s.
- The port's depends on it a little: 9.5 s for 310,759 annotated records.

Where the time goes (HG008 somatic, `-log_level DEBUG`, 2 threads, 3.2-3.6 s wall, 5.1 s CPU):
- the input VCF (50 MB bgzipped, 346,391 records, 17,254 kept): 0.45 s;
- then, in parallel:
  - the Ensembl cache, 0.7 s;
  - ClinVar (2.2 million records scanned), 1.2 s;
  - the mappability file (1.0 GB of text, 45 million lines), 2.4 s, which is the critical path;
- the chromosome tasks (impacts, phasing, annotation, encoding): 0.1 s;
- writing the BGZF file and index: 0.2 s.

With file-based PON and gnomAD (synthetic, 2 threads):
- the PON scan (764 MB, 25.7 million lines): 2.3 s;
- mappability: 2.5 s;
- gnomAD: 0.1-0.5 s per chromosome file;
- ClinVar: 1.1 s;
- processing 310,759 records: 2.9 s;
- writing: 2.0 s (single-threaded deflate at level 5).

**Scaling:** 5.4 s at 1 thread and 3.2 s at 2. From 4 threads, the mappability and PON files are inflated
on threads of their own beside the pool (gz `open_threaded`). That takes the mappability scan from 2.4 s
to 1.4 s, in one 4-thread run made to check that it gives the same output.

**Memory:** the scans keep only the entries at the input's positions. So memory stays at 0.2 GB on the
pipeline inputs, 0.9 GB with 310,759 records, and 1.1 GB when the transcript file's 677,000 lines are held
for writing.

## 5. Open points

- **GRCh37** (`-ref_genome_version 37`) is not tested; no 37 resources were at hand.
- **A plain VCF in Java mode** reproduces Java's failure (PA3), but not the contents of its truncated
  output.
- **Malformed resource lines:** the resource readers parse only what the lookups need. A malformed line
  elsewhere in the PON, gnomAD or mappability file passes in the port, where Java exits. ClinVar lines
  are parsed only at the positions needed.
- **A .vcf.gz input without an index:** Java exits ("An index is required"); the port reads it.
- **Sorted resource files:** the mappability and ClinVar lookups assume position-sorted resource files,
  as Java's searches do. The mappability reader warns if its file is not sorted.
- **Transcript file memory:** the transcript file's lines are held until the end (1.1 GB for 677,000
  lines). They could be written as each chromosome finishes.
- **Threads:** runs above 2 threads were not timed; the pipeline gives PAVE 8.

## 6. Shared code

PAVE added to `hmf-common-rs`:
- **`gz`:**
  - `GzReader`: zlib-rs inflate, one gzip member after another, moved from esvee-rs `caller::vcf`;
  - `LineReader`: lines as byte slices, ending in '\n' or "\r\n";
  - `open_threaded`: inflating on a thread of its own.
- **`vcf`:**
  - VCF headers: the standard line repair, the sorted output order, INFO and FORMAT definitions;
  - records: AbstractVCFCodec parsing and the getAttributeAs* accessors;
  - the VCFEncoder, from esvee-rs `caller::vcf`;
  - `VcfFileWriter`: the BGZF stream (65,498-byte blocks, zlib level 5) with htsjdk's on-the-fly tabix
    index, from amber-rs `vcf`. It includes htsjdk's check that positions do not go back within a
    sequence. PURPLE and LINX can write their VCFs with it.
- **`regions::positions_within`** (BaseRegion.positionsWithin).
- **Ensembl accessors:** mutable access to the Ensembl cache's gene and transcript maps, which LINX added
  too, for tests that fill a cache.
- **Dependencies:**
  - flate2 with the zlib backend, so that deflate gives Java's bytes;
  - memchr;
  - zlib-rs.

  Ports that already used flate2's default backend (SAGE for its gzip side files, REDUX for the jitter
  table) now get zlib's deflate through Cargo's feature unification. Their files' contents are
  unchanged.

## 7. Following upstream (883ab11bba, reviewed 2026-11-28)

The port started from ca04e6c39d. Upstream master 883ab11bba changes PAVE once; hmf-common's changes in that range
(PipelineToolDirectories, VirusInterpreterDataLoader, from f5f660cc9f) are code PAVE does not use.

| Upstream commit | Change | Decision | Reason | Evidence |
|---|---|---|---|---|
| 4adc014589 | Special splice regions: a built-in resource `splice/special_splice_regions.tsv` (gene, transcript, exon rank, exon start or end, bases from / to), read once; ImpactClassifier marks a variant overlapping one (on the listed transcript and exon, or the next exon) as in a splice region, as for the standard region; effects, HGVS and splice impact unchanged | accepted | The purpose is stated (AUS-424) and the one entry is right: ENST00000397752 is MET's canonical transcript (NM_000245, + strand) in the Ensembl cache, its exon 14 chr7:116,771,849-116,771,989 is the 141-base juxtamembrane exon (phases 1 to 1, an in-frame 47-codon loss) whose skipping is METex14; on the + strand the exon start is the acceptor, so bases -30 to -10 are intron 13's polypyrimidine tract and branch point region, where METex14 deletions and substitutions lie. The flag reaches the IMPACT field's canonical splice region (and with it PURPLE's driver candidates: SPLICE_REGION likelihood, not reported) and the transcript file's SpliceRegion column. It stays with that transcript and exon (MET's other transcript ENST00000318493, which has the same exon, is not flagged). Two edges are kept as written and asked about (PA10). | Java at ca04e6c39d and at 883ab11bba on the synthetic `met14` case: 137 of 323 records change, only the canonical splice region flag of ENST00000397752 (21 SNVs at -30 to -10, 3 MNVs, 113 INDELs); the port in Java mode is byte-identical to the new Java. No other case has a variant there: all 29 regression cases are unchanged. |
| 4adc014589 | SpliceContextTest.testSpecialSpliceRegion | accepted | | ported (`tests::splice_context::test_special_splice_region`) |

**References and results.** The 29 cases of `work/pave/scripts` and the new `met14` were run with the 883ab11bba jar
(`baseline/pave/java_883ab11bba/<case>`, same commands and threads as before) and with the port in Java mode
(`baseline/pave/rust_up/<case>`): every case that was identical before is identical (VCF, BGZF bytes, tabix index,
transcript file); the four that end with exit 1 in Java (plain_vcf_chr21: PA3; hg008_germline_reportable and
reannotate: PA9; pa4_complex: PA4) end so in the port with the same partial VCF, as before (the partial outputs
compare as they did: the VCF text equal, the truncated index files not). Times (2 threads, a loaded machine): HG008
somatic Java 20.8 s and 5.5 GB, the port 3.4 s and 0.20 GB; the synthetic set 87.6 s and 9.0 GB against 10.5 s and
0.90 GB; `met14` (1 thread, transcript file) 112.5 s and 8.6 GB against 6.7 s and 0.16 GB. Scripts: `work/up_pave_qsee/scripts/pave/` (cases.sh, run_one.sh, java_all.sh,
rust_all.sh, compare.sh, make_met14.py).

**Effect in the default mode** (the merged port, `rust_before_default`, against this one, `rust_up_default`): no
change on the 29 cases (none has a variant in MET's intron 13 region); on `met14` the 137 records above, as in Java.
