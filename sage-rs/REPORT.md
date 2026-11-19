# sage-rs: progress record

A running record of what was built, how it was checked, how it performs and what was found in the
Java code along the way. It is written so that it can later be turned into an announcement, an
upstream issue list or a pull request description. Detail lives in [PLAN.md](PLAN.md) (status and
measurements) and [QUIRKS.md](QUIRKS.md) (every Java behaviour reproduced on purpose).

Baseline: Java SAGE 5.1 at upstream commit `a224a2c`.

## 1. What exists

`sage-rs` is a Rust reimplementation of SAGE 5.1. It takes the same command-line options and writes
the same VCF, byte for byte, as Java run with `-threads 1`.

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

- **Definition.** The decompressed VCF must equal Java's with `diff`. Java is run with one thread
  because its output depends on thread scheduling (section 5); the port gives the same bytes at any
  thread count.
- **Data.** A private Illumina WGS tumor/normal pair (whole genome), and the public GIAB HG008
  tumor/normal pair for SBX and Ultima (chr21, after Java REDUX).
- **Whole genome.** 360,440 records without TINC and 611,876 / 365,286 with `-run_tinc`. Against Java
  at 16 threads the only differing lines are those with thread-order `LPS` ids (4,573 and 4,704).
- **Regression suite.** 30 comparisons against stored Java output, run after every change: each mode
  on each platform, MSI, multi-sample, append, CRAM. It is kept outside the repository because the
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
| With `-run_tinc` | 3,407 s | 764 s (first complete version), 593 s (second round), 327 s and 1.8 GB (current) |
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

- **mimalloc**: 10 to 15% faster, but memory grew through a whole genome (12 GB against 2 GB), because
  region results are allocated on worker threads and freed on the writer thread.
- **Slice comparison in the repeat search**: slower than the byte loop for the short repeats involved.

### Where the time goes now

Illumina chr21, one thread: BAM decompression and reading 34%, matching reads to candidates and
quality scores 13%, allocation 11%, decoding records 8%, jitter check 7%, read-context building 6%,
fragment sync 5%. For Ultima, read-context building is about half.

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

## 5. Issues found in the Java tools

These are candidates for upstream reports. The port reproduces every one that affects SAGE's output
(QUIRKS.md has the full list with file references). The ones that matter most to a user:

### Output depends on the thread count

- `LPS` ids come from a counter shared by all region threads, so a multi-threaded run does not
  reproduce itself. On chr21, 8 threads against 1 thread differs on 152 lines.
- `NEAR_HOTSPOT` is set after the region has been handed to the writer, so the flag can be missing.
- For contigs shorter than `threads` x 100 kb (chrM) the partition size depends on `-threads`, and
  the calls change with it.

### CRAM

- **CRAM and BAM of the same reads give different calls.** CRAM normally omits `NM`; htsjdk does not
  regenerate it; `NumberEvents.getOrCalcNm` then counts mismatches itself, and that count is not always
  the aligner's. Map-quality penalties change. Illumina chr21: 6,023 records from BAM, 5,908 from CRAM.
- htsjdk in SAGE 5.1 rejects CRAM 3.1, which current samtools writes by default.

### Crashes

- SAGE without `-sequencing_type` on the public SBX chr21 slice: null pointer in the candidate stage,
  after "ignoring read with position before prior position".
- REDUX 2.1: `BaseQualityResults.addPerfCounter` merges a counter into itself and fails during BQR.
- REDUX 2.1 with `-form_consensus` on the public Ultima slices: fails on a consensus read with an empty
  base.
- Ultima reads without `tp` / `t0` tags: null pointer rather than a message.
- Visualiser: `-vis_variants` with a variant two or more regions out of position order stops the VCF
  writer ("Features added out of order").
- Visualiser: with `-vis_purple_vcf`, a reportable variant whose transcript has no amino acid entry
  (intergenic, non-coding), or whose protein change is written in a form `AminoAcidEvent.parse` does
  not expect (`p.Gly12Ter`), is a null pointer exception: a warning, and no file for the variant.

### Logic that looks unintended and changes calls

- `PhasedGroupCollection.merge` uses `putAll`, which drops phase groups on a key collision.
- `PhasingUtils.mergeUninformative` tests a group against itself, so the test is always true.
- `VariantFilters.belowMinTumorVafProbability` rounds a ratio to 0 or 1 before multiplying.
- `FragmentSync`: when the merge throws, the second mate is dropped and the first processed alone.
- `tinc.VariantData.isMnv` requires alt to equal ref, so it is never true.
- `TincCalculator.calculate` has an early exit on a counter that is never incremented.
- `ReadContextCounterFactory` calls the getter `markedInvalid()` where `markInvalid()` was meant.
- `BedFileReader.loadBedFileChrMap` loses the earlier block of a BED file that revisits a chromosome.
- `IndelDeduper.checkDedupCombinations` orders tied combinations by identity hash code, so which of
  two variants is marked `dedupIndel` depends on JVM state. Observed: the same site gives opposite
  answers when the run's region starts 1.2 Mb earlier, or when a second reference sample is added.
- Several results depend on `HashMap` iteration order (leftover cached mates, uninformative phase
  sets, the Ultima maximum raw quality, read-length sampling).
- The binomial tail `1 - cdf` is exactly 0 in some cases only because of `FastMath` rounding, which
  selects a different penalty; one `MQF` value in a whole genome depends on it.

- **Turning the visualiser on removes every chrM call.** `VariantVis` needs the chromosome length from
  `RefGenomeCoordinates`, which lists 1-22, X and Y; on chrM the counters cannot be built and the VCF
  has no chrM records (0 against 128 on the Illumina pair), with nothing in the log.
- With `-vis_purple_vcf` and no reportable variant, SAGE exits leaving an empty, unterminated
  `.vcf.gz`.

### Visualiser output

- The reads shown, and the order of reads that tie, come from `ThreadLocalRandom`: two runs never
  give the same file.
- Text positions come from AWT font metrics, so the files depend on the fonts of the machine.
- Several filters are listed in an order that partly depends on identity hash codes.
- `BaseSeqViewModel.create` indexes the mates of a merged fragment by read-base index where a
  reference position is meant, so the overlap border is misplaced after an indel in the merged read.
- The `AED` row shows the first sample's two averages (all reads, alt reads) under "Tumor" and "Normal".
- `MQF` is rounded to a whole number (`round(round(x * 10) / 10)`).
- After a frameshift reaches a stop codon, `AminoAcidUtil` uses an index into the fetched bases as a
  genome position for the rest of the track.
- An unknown three-letter amino acid code in a `delins` is drawn as the letters of `null`.

### Cosmetic

- `RC_JIT` is declared with three values and written with two.
- The `fragmentLength` filter is declared but its test is commented out.

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

## 7. Not done

- Split reads compared with Java on real RNA data.
- Ultima reads formed by REDUX consensus (REDUX fails on the test data, above).
- Further speed: decompression (34 to 40% of CPU) and allocation (about 10%) are what is left; the
  visualiser builds its files on the writer thread.
