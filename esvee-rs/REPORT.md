# ESVEE in Rust: record of the port

A record kept for a later announcement and for reports to the Java maintainers. Status is in
[PLAN.md](PLAN.md); how to build and run is in [README.md](README.md).

Java reference: ESVEE 2.1 built from upstream master (commit `ca04e6c39d`; a fix for prep BAM loading
that is needed is in it). Paths of Java files below are under
`esvee/src/main/java/com/hartwig/hmftools/esvee/` unless they start with `hmf-common`. Folder names such
as `baseline/esvee/...` and `tools/esvee/...` are relative to the root of the workspace in which the
port was made; they are outside this repository, because the outputs are derived from patient data.

## 1. Summary

- `esvee-rs` holds the four steps of ESVEE 2.1: `prep`, `assembly`, `depth`, `caller`. Each takes the
  Java step's command line and writes the same files.
- On a small Illumina tumor/normal input, every output of the Rust chain run at 8 threads equals the
  output of the Java chain run as a deterministic reference (section 2.5): both prep BAMs and the three
  prep TSVs, `assembly.tsv`, `breakend.tsv`, and the raw, ref-depth, unfiltered, somatic and germline
  VCFs.
- Each step was also compared with Java on its own, over option cases and other kinds of input: 42
  cases for prep, 52 for assembly, 93 paired runs for depth, 96 cases for the caller. Public GIAB HG008
  SBX and Ultima reads were used for prep and assembly, and an exome pair for prep and assembly.
- A deterministic Java reference had to be defined first. Java's prep output depends on the thread
  count, and Java's assembly output changes from run to run even at one thread (section 3). The
  reference is Java at one thread with every identity hash made equal (`-XX:hashCode=2`). The port
  gives that result at any thread count.
- 83 findings about the Java code are recorded in section 4: 11 in prep, 36 in assembly, 19 in depth,
  17 in the caller. The port reproduces those that change the output.
- On the small input at 8 threads the Rust chain takes 24 s against Java's 95 s, with a peak of 1.9 GB
  against 44 GB (section 5).
- On an Illumina whole-genome tumor/normal pair (4.2 billion reads) every output of the Rust chain
  equals the Java reference (section 2.6). The prep tables are identical, and the prep BAMs hold the
  same records; in each BAM one pair of tied neighbouring records is in the other order, which is
  Java's own nondeterminism (finding P7). From those prep outputs on, `assembly.tsv`, `breakend.tsv`
  and all five VCFs are identical, whether the Rust steps read the Rust or the Java prep output.
- On that pair the Rust chain at 8 threads takes 14 minutes against 2.9 hours for the Java reference
  at one thread. The assembly step peaked at 33.4 GB at first; it now peaks at 7.2 GB with the same
  output (section 5).
- Not ported: htsjdk's record validation behind `-bam_validation` in the
  depth step. Not tested: reference genome 37 on real data, in any step (section 6).

## 2. Verification

What "identical" means throughout: TSVs byte for byte; VCFs byte for byte after decompression, header
included; BAMs as the same records in the same order with the same fields and attributes, header equal
apart from `@PG` lines. Compressed bytes and index files are not compared; indexes were checked by
querying through them.

Inputs:
- **Small pair:** slices of an Illumina whole-genome tumor/normal pair after REDUX, 120 million reads in
  the two BAMs; and a smaller cut of it for option checks.
- **Exome pair:** an Illumina exome tumor/normal pair of another sequencing run, all chromosomes, 0.9
  billion reads.
- **Public SBX and Ultima:** chr21 of GIAB HG008 tumor and normal after REDUX, as BAM and as CRAM 3.0.
- **Synthetic:** generated inputs holding nothing from real data (depth and caller steps; SAGA
  resources for prep and assembly).

### 2.1 Prep

Java at one thread against the port. Scripts: `tools/esvee/prep_cmp.sh`, `prep_optcmp_all.sh`.

| Input | Runs | Result |
|---|---|---|
| Small pair | Rust at 8 threads, at 1 thread, with sort runs of 60 MB, with debug logging | three TSVs identical; 770,828 + 627,454 BAM records in the same order with the same attributes |
| Smaller cut of the small pair, Rust at 8 threads | default; tumor only; normal only; tumor only with `-no_trim_read_id`; `-no_cache_bam` (pair and tumor only); `-specific_chr` (two and one); `-specific_regions` (both syntaxes); `-write_types` with PREP_READ (`read.tsv` identical), PREP_STANDARD;PREP_READ with `-no_cache_bam`, PREP_JUNCTION only, PREP_BAM only, PREP_READ only; `-partition_size` 500000 and 3000000; `-min_align_bases 60 -min_junction_frags 3`; `-min_junction_frags 1`; `-max_frag_length_override 600`; no `-bamtool` (unsorted BAMs in the same order); `-output_id`; no hotspot file; `-sequencing_type SBX` and `ULTIMA` on Illumina reads; `-no_clean_up`; `-perf_debug`; `-log_level DEBUG -bam_validation SILENT`; `-ref_genome_version 37` (no reads match: empty outputs); `-log_read_ids`; `-saga_fasta` with a synthetic resource (pair, and tumor only with `-no_cache_bam`); Rust at 1 and 3 threads | all identical; with `-perf_debug -log_debug` Java stops with an exception (finding P9) |
| Public SBX pair | `-sequencing_type SBX`; with `-specific_chr chr21`; tumor only; tumor only with `-no_cache_bam` | identical |
| Public Ultima tumor | `-sequencing_type ULTIMA`; with `-specific_chr chr21` | identical |
| CRAM 3.0 | public SBX normal; Illumina chr21 pair with `-specific_chr chr21`; Illumina tumor alone | identical |
| Exome pair, Rust at 8 threads | default | three TSVs identical; the same 780,908 + 544,687 BAM records; 3 pairs of neighbouring records in the tumor BAM and 2 in the normal BAM are in the opposite order |

- The final pass with the finished binary against the stored Java outputs: 41 cases identical and the
  one case in which Java exits with an error.
- The swapped pairs of the exome run are Java's own nondeterminism (finding P7). Each pair has the same
  position and strand, and both reads were claimed from the cache BAM. A second Java one-thread run on
  the same input differs from the first in the same way: the normal BAM identical, the tumor BAM with
  2 such pairs swapped.
- Rust at 8 threads and at 1 thread write byte-identical BAM files, as do two 8-thread runs.
- CRAM input first differed because htslib adds `MD` and `NM` attributes that the CRAM does not store
  and htsjdk does not add; the port turns that off.
- Unit tests: 31, of which 21 are the Java tests of the package; the rest check the emulated HashMap
  order against a JVM and the record helpers.

The whole genome is in section 2.6. Not tested: sambamba as `-bamtool` (treated like samtools;
not installed, so not compared); `ReadIdTrimmer` switching itself off part-way through a partition
(ported as written, no test data reaches it); reference genome 37 with reads that match it.

### 2.2 Assembly

Java at one thread with `-XX:hashCode=2` against the port at 8 threads. Script:
`tools/esvee/assembly_optcmp.sh all`. Each case compares every output file and the log text (lines that
are not DEBUG, sorted). With `-write_types ASSEMBLY_ALL` the files are `alignment.tsv`, `assembly.tsv`,
`assembly_read.tsv`, `breakend.tsv`, `decoy_match.tsv`, `phase_group_building.tsv`,
`phased_assembly.tsv` and the raw VCF.

All 52 cases are identical (41 without the visualiser, 11 with it):

| Kind | Cases |
|---|---|
| Small pair, modes | default; `-write_types ASSEMBLY_ALL`; tumor only (with all types and with the default types); germline only (reference sample alone, on the output of a normal-only prep run); `-esvee_prep_dir`; Rust at 3 threads |
| Small pair, options | `-specific_chr`; `-specific_regions`; `-specific_junctions` (also with `-run_assembly_vis`); `-purple_dir`; `-write_candidate_reads`; `-assembly_detailed_tsv`; `-disc_only_disabled`; `-asm_map_qual_threshold` 0 and 30; `-phase_process_limit`; `-output_id`; `-remote_phase_read_check_threshold`; `-disc_rate_increment`; `-asm_ref_base_write_max`; `-write_types` VCF, BREAKEND;ALIGNMENT, JUNC_ASSEMBLY;PHASED_ASSEMBLY;ASSEMBLY_READ and NONE; `-log_debug`; `-perf_debug -perf_log_time`; `-asm_build_debug -asm_log_all_reads`; `-log_read_ids` |
| Sequencing type | public SBX tumor and normal with `-sequencing_type SBX`; public Ultima tumor alone with `-sequencing_type ULTIMA`; the Illumina input run as SBX and as Ultima; the Illumina input prepared and assembled as SBX |
| Exome | the exome pair |
| SAGA | a synthetic resource built from the deletion and insertion breakends of the reference run, 381 variants including near misses (523 assemblies and 554 breakends matched); tumor only with `-assembly_detailed_tsv`; prep output carrying SAGA matches; an invalid resource (both tools stop the same way) |
| Decoy | `-decoy_genome` with a BWA image |
| Visualiser (`-run_assembly_vis`) | junctions of 2 and of 8 breakends of each SV type (10 and 39 files); assemblies with several breakend pairs (4 files); `-output_id` (10); `-purple_dir` (10, the same files as without it); `-log_debug` (5, DEBUG lines of the visualiser identical); public SBX pair (18); tumor only, germline only and public Ultima tumor (Java writes no files, findings A31 and A32); indel-only assemblies (no breakends, no files); the flag without `-specific_junctions` (ignored) |

- **BWA.** `alignment.tsv` holds every raw field of each BWA record (location, coordinates, MAPQ,
  cigar, score, flags, NM, XA, MD). It is identical to Java's in every case, about 22,000 alignment
  records over the Illumina, SBX, Ultima, exome, germline and SAGA runs. The decoy and SAGA results,
  which use two other indexes and parameter sets, are identical too.
- The port at 1 thread and at 8 threads gives the same files.
- The public SBX case found one difference during the work: Java's `DoubleStream.sum()` is a
  compensated sum (finding A29). The port now sums the same way.
- Unit tests: 101 ported from all 20 Java test classes of the package (108 with those of the shared
  modules). The emulated identity-hash set is checked against real Java output for 60 random sequences
  of adds and removes.

Not tested: reference genome 37 (the only genome 37 prep output at hand has no junctions); a real
decoy genome image (none was available: the SAGA test image stood in for it).

**Visualiser.** Each case compares the `vis/` folder by file name and each file byte for byte: 96
files in the 7 cases where Java writes any (some assemblies recur between cases). Two Java runs of one case wrote identical files, so there is
nothing random to allow for. Text positions come from AWT font metrics, so Java's files depend on the
machine's fonts; the port holds the metrics of the machine it was compared on (macOS). The drawing
code (`src/common/vis`: SVG as JFreeSVG writes it, the HTML that j2html renders, the hmf-common view
models) is checked by 15 unit tests against values dumped from the Java classes on synthetic input.
The Java package has no unit tests. The port takes 0.5 to 5 s on these cases where Java takes 65 s to 34
minutes; a run without the flag is as fast as before (2.7 to 2.9 s on the small pair).

### 2.3 Depth

Java against the port, both at the same thread count. Scripts: `tools/esvee/depth_checks.sh` (the list
of every comparison), `depth_regress.sh`.

93 paired runs: 85 identical, 8 with documented expected differences, none unexpected. 68 of the runs
write a VCF in both tools, and all 68 VCFs are identical, each with a tabix index that answers a query:
46 on the real input and 22 on synthetic inputs.

Real input (the raw VCF of the small pair, 2,425 variants, and its two BAMs):
- threads 8, 4, 3, 1, 0 and the option absent;
- a second input VCF (Java assembly at one thread, which gives a different VCF), at 8 and 1 threads;
- `-vaf_cap` 0.05, 0.5, 0, -1; `-proximity_distance` 0, 500, 100000;
- `-input_vcf` absent (found from the other options), `-output_id`, `-output_dir` absent;
- `-unmap_regions` absent; `-specific_regions` (one, three in both syntaxes, invalid, matching
  nothing); `-specific_chr`;
- samples in swapped order, tumor only, normal only, tumor only with `-unmap_regions` (both tools fail);
- `-ref_genome_version` 38, V38, HG38, 37 (no chromosome matches: header-only output), 39 (error);
- `-vcf_tag_prefix` on the raw VCF (both fail), on a VCF that already has the depth tags, with
  `-vaf_cap`; re-annotating a depth VCF without a prefix;
- `-bam_validation` STRICT, LENIENT, SILENT (no effect on these BAMs in Java), and an invalid value;
- `-log_debug`, `-log_level` DEBUG, warn, ERROR and invalid, `-perf_log_time`;
- CRAM 3.0 of one chromosome at 8 and 1 threads; CRAM 3.1 (Java cannot read it);
- errors: each required option missing, missing BAM, missing VCF, unknown option, stray argument, a
  value starting with `-`, non-numeric values, one BAM for two samples, a sample not in the VCF,
  `-help`, `-version`.

Synthetic inputs (random reference; reads with soft clips, deletions, insertions, skips, zero reference
length, duplicates, secondaries, supplementaries with `SA` in several forms, unpaired, mate unmapped,
mate elsewhere, three records of one name; dense variant clusters; small variant fractions so that caps
are hit):
- three standard sets, each with four option combinations: identical. Some contigs there have a
  sequence index above 127 (finding D1);
- the same reads and variants with those contigs at an index below 128: identical to Java, and Java's
  own output differs from the first form (finding D1);
- a chromosome in two blocks, records on a non-chromosome contig, genotypes with few or no values:
  identical; with unmap regions both tools fail (finding D9);
- unusual VCF text (header lines needing repair, empty and flag INFO, repeated filters, GQ/PL/FT,
  phased GT, several ALTs, missing ID, QUAL forms): identical;
- variants at the excluded poly-G positions on every chromosome: identical (finding D3);
- an unsorted BAM and a malformed `SA` attribute: both tools stop with status 1;
- `-perf_log_time` with every slice over the threshold, `-log_level TRACE`, `-bam_validation LENIENT`:
  identical.

The expected differences are of these kinds: both tools stop with status 1 and only the wording of the
message, or what is left on disk, differs (Java leaves a truncated VCF where the port leaves a complete
one); `-bam_validation STRICT` on a BAM htsjdk finds invalid, which is not ported; TRACE "ref limit ...
reached" lines (finding D19); CRAM 3.1, which Java cannot read.

Unit tests: 20, among them the two tests of the Java `DepthAnnotationTest` with the same cases and
assertions.

The whole genome (52,061 variants) is in section 2.6. Not tested: a real GRCh37 BAM
(`-ref_genome_version 37` was run on genome 38 data, where no chromosome name matches, and on the
synthetic poly-G set). Malformed input on which Java throws inside htsjdk is not matched case by case
(for example a mapped read without a cigar, or an `SA` attribute that is not a string): the port either
carries on or stops with its own message.

### 2.4 Caller

Java against the port. Scripts: `tools/esvee/caller_optcmp.sh`, `caller_synth.sh`. Each case compares
every VCF, `final_breakend.tsv` where written (the filter column as a set, finding C1), and the exit
code.

- **Standard run** on the ref-depth VCF of the small pair: unfiltered, somatic and germline VCFs
  identical (2,425 / 31 / 1,515 records).
- **Option matrix, 96 cases, all pass:** tumor and reference, tumor only, germline only;
  `-write_breakend_tsv` in each mode; `-output_id`; default input and prep paths taken from
  `-output_dir`; artefact PON files (both, SV only, SGL only, without a main PON); `-pon_margin` and
  `-pon_sgl_margin` (larger, zero, negative); PON in TSV format; another PON; no PON; an unordered PON
  (exit 1); `-target_regions_bed` (with `-min_qual`, tumor only with `-filter_sgls`, with SBX, a
  shuffled BED, a BED with a header and gzip); `-min_qual`, `-hotspot_min_qual`, `-min_support`,
  `-hotspot_min_support`, `-sgl_min_support`, `-filter_sgls`; `-sequencing_type` ILLUMINA, SBX, ULTIMA
  and invalid; `-manual_ref_depth` on the VCF before depth annotation (paired, tumor only, germline
  only) and on one that has depth; `-specific_chr`; `-specific_regions`; without repeat mask, hotspots
  or any resource; an uncompressed repeat mask; `-ref_genome_version` 37, V38, invalid and missing;
  `-log_debug`, `-log_level`; an uncompressed input VCF, a header-only VCF, a breakend without its
  mate; a known-pair file built from the input so that 547 breakends are hotspots (with the hotspot
  thresholds, target regions, `-filter_sgls`, SBX); sample names that are wrong, swapped, identical,
  missing, or only a substring of the VCF's; missing prep files or VCF; malformed PON, hotspot and
  repeat-mask lines; an unknown option, an option without a value, a non-numeric value.
- Without `-output_dir` (files prefixed `null` in the working directory, finding C8): run by hand,
  identical.
- **Synthetic data** with made-up coordinates: identical to Java in five modes; the same comparison
  runs in `cargo test` against the stored Java output.
- **Unit tests:** 38 pass and 1 is ignored, as it is in Java. All 18 Java tests of the package are
  there under their names. They were ported from the source; they were not run in Java here.

Not covered by the real test data, only by unit tests or the synthetic set: the filters `strandBias`,
`lineSource`, `sbxInvZeroMotif` and `invShortFragLowVaf`; SAGA-matched variants; contigs outside
chr1-22/X/Y. `-sequencing_type SBX` and `ULTIMA` were run on Illumina data only. A reference genome 37
input was not available; only the case of version 37 given for genome 38 data was run.

### 2.5 End to end

The four steps were run in order on the small pair, each step reading the previous step's output of its
own chain:

- Java reference chain: every step at one thread (`tools/esvee/run_java_ref.sh`, output in
  `baseline/esvee/java/small_ref`). The prep output is that of an earlier Java one-thread run; assembly,
  depth and caller were run from it with `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`;
- Rust chain at 8 threads (`tools/esvee/run_rust_steps.sh`, output in `baseline/esvee/rust/small_t8`);
- comparison of every output: `tools/esvee/cmp_all.sh`.

| Output | Result |
|---|---|
| `prep.junction.tsv`, `prep.fragment_length.tsv`, `prep.disc_stats.tsv` | identical |
| normal prep BAM | same records, same order (770,828) |
| tumor prep BAM | same records, same order (627,454) |
| `assembly.tsv`, `breakend.tsv` | identical |
| raw VCF | identical (2,425 records) |
| ref-depth VCF | identical (2,425 records) |
| unfiltered VCF | identical (2,425 records) |
| somatic VCF | identical (31 records) |
| germline VCF | identical (1,516 records) |

The germline VCF has one record more than in section 2.4 because the input chain differs: there the
caller read the output of a Java chain run at 8 threads.

### 2.6 Whole genome

An Illumina whole-genome tumor/normal pair after REDUX, 4.2 billion reads in the two BAMs (tumor 120 GB,
normal 96 GB).

- Java reference: prep at one thread; assembly, depth and caller at one thread with
  `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, each reading the previous Java step's output
  (`baseline/esvee/java/wgs_t1`).
- Rust chain at 8 threads (`baseline/esvee/rust/wgs_t8`).
- Rust assembly, depth and caller at 8 threads on the Java prep output
  (`PREP_IN=... run_rust_steps.sh ... assembly`, output in `baseline/esvee/rust/wgs_javaprep_t8`), to
  compare the later steps on the same input.

| Output | Rust chain | Rust on Java's prep output |
|---|---|---|
| `prep.junction.tsv`, `prep.fragment_length.tsv`, `prep.disc_stats.tsv` | identical | (Java's) |
| normal prep BAM | same 20,377,496 records; one pair of neighbouring records in the other order | (Java's) |
| tumor prep BAM | same 16,658,562 records; one pair of neighbouring records in the other order | (Java's) |
| `assembly.tsv` (85,306 rows), `breakend.tsv` | identical | identical |
| raw, ref-depth, unfiltered VCF | identical (52,061 records each) | identical |
| somatic VCF | identical (3,072 records) | identical |
| germline VCF | identical (35,338 records) | identical |

- Each swapped pair has the same position, strand, flags and alignment, so the two orders are equally
  sorted. Their read names agree in the lowest 16 (normal) and 14 (tumor) bits of Java's spread string
  hash, which puts them in the same bucket of a hash table; two random names agree that far about
  once in 10,000 pairs. That is finding P7: Java orders such reads through a set keyed by identity
  hash, and two Java runs on the exome pair already differed in this way. The swap changes nothing
  downstream: the Rust chain's later outputs equal Java's.
- The same Java chain run first without `-XX:hashCode=2` (default identity hashes; kept in
  `wgs_t1/plain_hash`) differs from the reference in every output from assembly on. Counting a record
  by chromosome, position, ref, alt and filter: 14 rows of `assembly.tsv` differ; the raw VCF has
  52,062 variants instead of 52,061 (3 not in the reference, 2 of the reference missing); the somatic
  VCF 3,078 calls instead of 3,072 (6 more, none missing); the germline VCF 35,337 instead of 35,338
  (11 not in the reference, 12 missing). Java's somatic call set thus depends on the identity hashes,
  which differ in every JVM run (section 3). Only one such run was made on this input.

## 3. Determinism of Java

A comparison needs a Java output that is the same on every run. Measured on the small pair:

| Step | One thread, twice | 8 threads, twice | 8 threads against one |
|---|---|---|---|
| prep | all outputs identical (but see P7 for the exome pair) | differ | differ |
| assembly | differ | differ | differ |
| depth | | | VCF and index identical byte for byte |
| caller (no thread option) | the three VCFs and their indexes identical byte for byte over 10 runs; the `Filter` column of `final_breakend.tsv` differs between all four runs that wrote it | | |

### Prep depends on the thread count

- Two 8-thread runs: `junction.tsv` has the same rows in a different order (2,828 lines moved); the
  `PrepReads` count in `disc_stats.tsv` is 1,393,940 against 1,393,945; the BAMs differ by 1 record one
  way and 4 to 5 the other.
- 8 threads against one: `junction.tsv` same rows, different order (16,540 lines); `PrepReads` 1,393,940
  against 1,394,477; the 8-thread BAMs lack 181 (normal) and 184 (tumor) records that one thread
  writes, and hold 2 and 5 that one thread lacks. `fragment_length.tsv` is identical.
- Cause: the partitions share a cache of spanning reads, and a partition's lookup and update of it are
  two separate synchronised calls (finding P1). Junction rows are appended as partitions finish.
- At one thread the final outputs were identical between two runs of the small pair. The temporary cache
  BAMs (`-no_clean_up`) differ in record order even then, and on the exome pair so does the order of a
  few tied records in the final BAM (finding P7).

### Assembly changes from run to run

- Two one-thread runs: `assembly.tsv` identical. In `breakend.tsv` the `AssemblyId` differs on nearly
  all rows, `AvgFragmentLength` on 42, `IncompleteFragments` on 36, `Id`/`MateId` on 4. In the VCF
  `ASMID` and `SEGID` differ on nearly all records, `AVGLEN` on some, and 4 records differ in other
  fields.
- Two 8-thread runs: the same kind of difference, and also `FacingBreakendIds` on 2 rows and ids on 6.
- One thread against 8: in addition `PhaseGroupId` in `assembly.tsv` (4 to 9 rows), `DiscFrags` and
  `RefDiscFrags` (1 row between two of the 8-thread runs), and a log line (finding A3).
- Cause of the one-thread variation: iteration over `HashSet` and `HashMap` keyed by objects without
  `hashCode()` (`PhaseGroup`, `JunctionAssembly`, `JunctionGroup` and others), so the order is that of
  the identity hashes, which differ in every JVM run (finding A1).
- Proof: with `-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, which makes every identity hash 1, two
  one-thread runs are byte-identical in every output, with all write types, the compressed bytes of the
  VCF included.
- On the whole genome the identity hashes change the calls, not only ids: a one-thread Java chain with
  default hashes against one with `-XX:hashCode=2`, from the same prep output, has 6 more somatic
  calls (3,078 against 3,072) and 23 germline records not shared (section 2.6).
- With constant hashes, one thread against 8 still differs: 7 lines of `assembly.tsv`, 3 of
  `breakend.tsv`, 3 of the VCF, 339 of `assembly_read.tsv`; `alignment.tsv` and `phased_assembly.tsv`
  in row order only; `phase_group_building.tsv` in order and in 27 lines. One cause found: support
  reads shared between assembly alignments are read and overwritten by each (finding A2). Others lie
  in phase group building, where the result depends on the order in which groups are linked.

### The reference, and how the port reproduces it

The reference is Java with `-threads 1` in every step and `-XX:+UnlockExperimentalVMOptions
-XX:hashCode=2`. The option matters for the assembly step only.

- **Prep.** A partition is read and its junctions found on any thread. The stage that reads the
  spanning cache, updates it and writes output runs in partition order under one lock, so the result is
  that of one thread taking the partitions in order. Sets that Java keys by identity hash are kept in
  insertion order.
- **Assembly.** Junction groups, phase groups and assembly alignments are processed in parallel.
  Everything whose result depends on processing order is done in Java's one-thread order: phase group
  building, the allocation of support reads to breakends, decoy rows, TSV row order. Where Java
  iterates an identity-keyed set, the port uses a set that iterates as Java's `HashMap` does when all
  hashes are equal (one bucket, as a list and then as a red-black tree).
- **Depth.** Java is deterministic. The port hands out groups of variants rather than chromosomes to
  its threads; groups share no state in Java either, so the result is the same.
- **Caller.** The VCFs are deterministic. For `final_breakend.tsv` the port writes the filters in enum
  order and the comparison treats the column as a set.

## 4. Java findings

Numbered per step: P (prep), A (assembly), D (depth), C (caller). "Reproduced" means the port behaves
the same way so that its output equals Java's; where Java stops with an exception, "reproduced as a
failure" means the port stops with status 1 at the same point. Where no observed effect is given, the
finding comes from reading the code while porting and its effect on output was not measured; findings
the notes mark explicitly as not observed say so.

### 4.1 Prep

**Wrong result**

- **P3. `prep/DepthTracker.java:27`** `(int)ceil(regionLength / mWindowSize)` divides two integers, so
  `ceil` has nothing to round. The last partial window of a region (end of a chromosome, specific
  regions) has no depth, and junctions there get depth 0. Reproduced.
- **P4. `prep/SpanningReadCache.java:69`** names a partition `position / partitionSize`, while
  partitions start at `1 + i * size` (`prep/ChromosomeTask.java:131`). The last base of a partition is
  named as the next partition. Reproduced.
- **P5. `prep/CandidateBamWriter.java:195`** makes a new `PrepRead` for a read claimed from the cache
  BAM. Its filter flags are lost, so a read with a poly-G soft clip that
  `prep/ResultsWriter.java:256` (`filterBamRecord`) would hold back is written. Observed: Java's
  one-thread prep BAMs of the small pair hold 8 (tumor) and 6 (normal) such records; with
  `-no_cache_bam` they are not written. The two modes also claim candidates differently: the cache BAM
  mode registers a junction fragment's id against the chromosomes of its own remote partitions
  (`prep/CandidateBamWriter.java:54`) and claims a cached read only when it lies on one of them; the
  in-memory mode finds candidates through the partitions the candidate itself points to. Observed on
  the smaller cut: `-no_cache_bam` writes one supplementary read more (85,731 against 85,730 tumor
  records). Reproduced.
- **P6. `prep/CandidateBamWriter.java:77`** files an unmapped candidate read under chromosome "-1";
  line 192 and the loop over `mChrJunctionReadIds` never open that file, so such reads are never
  claimed. Reproduced.
- **P8. `prep/DiscordantGroups.java:263`** adds a remote discordant junction to the list of the
  partition that found it, and it is written with that partition's chromosome
  (`prep/PartitionSlicer.java:255`). Observed on the small pair: 8 such junctions lie outside their
  partition, on the same chromosome, so they are only out of order. A remote region on another
  chromosome, which hotspot pairs allow, would be written under the wrong chromosome: from reading the
  code, not observed. Reproduced.
- **P10. `prep/ReadFilters.java:107-113`** A supplementary read with fewer than 40 aligned bases is
  always dropped: the score function returns -1 for it, and the last line tests
  `-1 < MIN_SUPP_CALC_ALIGNMENT_SCORE`. The test at line 104 (`alignedLength < 40 && hasLowMapQual`)
  suggests that the map quality was meant to matter. Reproduced.

**Nondeterminism**

- **P1. Output depends on the thread count.** `prep/PartitionSlicer.java:121` reads the expected read
  ids before the partition's own cache update at line 242. Both are synchronised calls on the shared
  `SpanningReadCache` (`prep/SpanningReadCache.java:80`, `:87`), but not one step. A partition B whose
  candidate read's fragment has its junction in partition A loses the read when A finishes between B's
  two calls: A registered the id after B looked, and B then drops the candidate because A is already
  processed (`prep/SpanningReadCache.java:192`). `ResultsWriter.writeJunctionData`
  (`prep/ResultsWriter.java:156`) appends rows in completion order. Observed effect: section 3. Not
  reproduced: the port gives the one-thread result at any thread count.
- **P7. Identity-hashed set.** `prep/JunctionTracker.java:124` is a `Set<ReadGroup>`, and `ReadGroup`
  has no `hashCode`; it is iterated at line 188. The order of candidate groups within one hash bucket
  of the spanning map therefore varies between runs. Observed: the cache BAMs of two one-thread runs
  hold the same records in a different order; and in final output, two one-thread runs on the exome
  pair give tumor prep BAMs that differ in the order of 2 pairs of records sharing position and strand.
  Not seen on the small pair; on the whole genome one such pair in each prep BAM differs from the Java
  run (section 2.6). Not reproduced: the port keeps insertion order, and so can differ from a given
  Java run in the order of such pairs (sections 2.1 and 2.6). It did not change any later output on
  the whole genome.

**Crash**

- **P9. `-perf_debug` with debug logging stops Java.** `hmf-common` `PerformanceCounter.java:245`
  compares times with `first.Time < second.Time ? 1 : -1`, never 0. `Collections.sort` at line 172
  throws "Comparison method violates its general contract" (reached from
  `prep/ChromosomeTask.java:96`), exit code 1. Observed on the smaller cut. Not reproduced: the port
  completes.

**Oddity**

- **P2. `PrepReads` in `disc_stats.tsv`** (`prep/PrepApplication.java:76`) is the count taken before
  the cached candidate reads are written (line 81), so it is lower than the "records written to BAM"
  log line. Reproduced.
- **P11. Output depends on read names.** `prep/JunctionTracker.java:259` iterates a
  `HashMap<String,..>` keyed by read id. Renaming the reads can change which read is a junction's
  `InitialReadId` and the order of tied BAM records. Not a fault as such. Reproduced, by emulating
  `HashMap` order.

### 4.2 Assembly

**Wrong result**

- **A4. `assembly/phase/RemoteGroupBuilder.java:112`** A `Set<JunctionGroup>` uses
  `BaseRegion.equals/hashCode`, which look at start and end only. Junction groups with equal bounds on
  different chromosomes collapse into one, and the iteration order follows that hash. Reproduced.
- **A6. `assembly/types/Junction.java:223`** tests `otherSupportFragsIndex != null` but reads
  `values[otherJuncFragsIndex]`: the null test is on the wrong column index. Reproduced.
- **A9. `assembly/types/RemoteRegion.java:129`** `addSoftClipMapQual(1, nextRegion.softClipMapQualTotal())`
  is inside the loop over the three read types, so the merged region gets the other region's soft-clip
  map-quality total three times. Reproduced.
- **A10. `assembly/SeqTechUtils.java:103`** In the second branch of `trimSbxUncertainBases`,
  `trimCountStart` is computed with `fromStart` false. Reproduced.
- **A11. `assembly/SeqTechUtils.java:39`** `PROXIMATE_JUNCTION_DISTANCE` is set to the value for
  unpaired reads (100) for every sequencing type; the assignment stands outside the branches for the
  types. Reproduced.
- **A12. `assembly/JunctionGroupAssembler.java:286`** `trimIlluminaAdapterBases` is called before mate
  links exist (`mateRead()` is null), so the adapter trim never happens in the application, only in
  the unit test. Reproduced.
- **A13. `assembly/IndelBuilder.java:56`** The second soft-clip test repeats `junction.isForward()`
  where reverse was meant, so left clips at reverse junctions are never counted. Reproduced.
- **A14. `assembly/LineUtils.java:348,350`** `maxPostPolyAtLength` is assigned for each read, and the
  `max` at line 352 is taken against a variable that is always 0. The last read wins instead of the
  maximum. Reproduced.
- **A17. `common/CommonUtils.java:154`** (through `prep/ReadFilters.java:116`, `filterLowQualRead`):
  bases that are not ACGT are rewritten to N only up to the point where the low-quality count exits
  early. Reproduced, in prep and in assembly.
- **A18. `assembly/alignment/Breakend.java:328`** `compareTo` passes `otherBreakend().Chromosome`
  twice, so it never compares with the other breakend's mate chromosome. Reproduced.
- **A19. `assembly/phase/ExtensionCandidate.java:94`** `secondJuncOffset` is computed from `first`, so
  the two offsets always tie and the rule "favour a more precise link" never applies. Reproduced.
- **A20. `assembly/phase/PhaseSetBuilder.java:617-618`** `if(unmappedReads.isEmpty()) break;` ends
  `findUnmappedExtensions` at the first assembly without unmapped reads; `continue` may have been
  meant. Reproduced.
- **A22. `assembly/output/AssemblyWriterUtils.java:199`** compares `repeat.totalLength()` with
  `longest` (the repeat with the highest count) instead of `longestSequence`. Reproduced.
- **A23. `assembly/phase/UnmappedBaseExtender.java:353`** `baseCounts[consensusNucIndex] = i` stores
  the index of the read, not a count. Reproduced.
- **A26. `assembly/alignment/AssemblyAligner.java:499-520`** `hasSuspectExtension` loops over the
  assemblies but returns in the first iteration, so only the first assembly is examined. Reproduced.
- **A27. `assembly/alignment/AssemblyAlignment.java:371`**
  `currentSeqLength = nextAssembly.extensionLength()` is an assignment where the same variable is
  accumulated with `+=` elsewhere in the method. Reproduced.

**Nondeterminism**

- **A1. `assembly/AssemblyApplication.java:381`** `gatherAssemblies` iterates a `HashSet<PhaseGroup>`
  whose order is that of identity hashes. The ids of assembly alignments and the order of derived
  assemblies change from run to run at one thread. Observed: `AssemblyId` in `breakend.tsv`, `ASMID`,
  `SEGID` and `AVGLEN` in the VCF, and the other differences listed in section 3. The same class of
  problem exists wherever an identity-keyed `HashSet` or `HashMap` is iterated (`PhaseSetBuilder`, the
  `Set<Breakend>` of `AlignmentFragments`, and others). Reproduced for the JVM option `-XX:hashCode=2`.
- **A2. `assembly/alignment/AlignmentFragments.java:76,87-88,250`** Support reads shared by several
  assembly alignments (phase sets with secondary links to the same assembly) are read and overwritten
  by each. With more than one thread the read-level `InferredFragLength` and `BreakendSupport` in
  `assembly_read.tsv` depend on completion order. Observed: 344 rows on the small pair. Not reproduced:
  the port allocates support in assembly alignment order, which is Java's one-thread result.
- **A3. `assembly/phase/RemoteReadExtractor.java:59-62` with `assembly/phase/PhaseSetBuilder.java:119`**
  `resetCounts()` zeroes the "total" counters for every phase group, so the log line "remote reads
  extracted(N)" is the count of each task's last group. Observed: 0 at one thread, 629 and 635 at 8
  threads. Reproduced (the one-thread value).
- **A5. `hmf-common` `ReadIdTrimmer`** (the static `READ_ID_TRIMMER`, `assembly/AssemblyConfig.java:157`)
  is shared between threads without synchronisation. The first read name that does not fit the first
  read's prefix switches trimming off from then on, so read ids depend on processing order and on the
  thread count. Reproduced as written.
- **A30. Decoy rows in completion order.** Decoy matches (`AlignmentChecker`, called at
  `assembly/JunctionGroupAssembler.java:219`) are written by each task as it goes, so `decoy_match.tsv`
  and the decoy rows of `assembly.tsv` are in task completion order with more than one thread. Not
  reproduced: the port writes them in junction group order, the one-thread order.
**Visualiser**

- **A31. Visualiser fails on single breakends.** `assembly/vis/AssemblyVisualiser.java:717` throws
  `RuntimeException("Got a single breakend")` for an assembly alignment with a single breakend.
  `AssemblyApplication.runVisualiser` (line 474) catches it, logs at DEBUG and prints the stack trace
  to stderr, so no file is written for such an assembly and the run reports success. Reproduced: the
  port skips the same assemblies and logs the DEBUG line.
- **A32. Visualiser writes nothing without a reference sample.** `AssemblyVisualiser.java:215` calls
  `mConfig.ReferenceIds.get(0)`, which throws `IndexOutOfBoundsException` in a tumor-only run, for
  every assembly; likewise `TumorIds.get(0)` in a germline-only run. The exception is caught as in
  A31, so `-run_assembly_vis` gives no files in those modes. Reproduced.
- **A33. `assembly/vis/AssemblyVisualiser.java:553-562,606`** (visualiser) Only the soft clip on the
  breakend's side is taken off the alignment's cigar, but the assembly bases passed with it are the
  aligned part alone. For an alignment clipped at both ends (the middle of a chain) the remaining clip
  reads past the end of the string: `StringIndexOutOfBoundsException` at
  `BaseSeqViewModel.java:183`, and no file for that assembly alignment. Seen for 4 of 8 chains of the
  small pair and 3 of the SBX case (one more SBX chain fails two lines earlier, at `:179`, on a matched
  base); this is the cause as read from the code and the stack trace.
  Reproduced (no file).
- **A34. `assembly/vis/ReadViewModel.java:401-402`** (visualiser) `handleHomology` is passed `readBases`
  for both the bases and the base qualities, so reads drawn across a homology are shaded by the
  character codes of their bases, which are all above the cut-off: no quality shading. Reproduced.
- **A35. Visualiser: objects changed while drawing.** `AssemblyVisualiser.java:568,575` reverses
  `AlignData.cigarElements()` (the alignment's own list, or a `subList` view of it) and
  `ReadViewModel.java:387,389-398` reverses `Read.cigarElements()` (when the read has no recorded
  mismatches) and the array of `Read.getBaseQuality()`, all in place. The step's outputs are already
  written, but an alignment that is the first segment of two breakends, or a read drawn twice, is
  seen changed the second time. The port keeps these lists per object for the length of the
  visualiser run, as Java does. In the 96 files compared it made no difference: a build that started
  from the unchanged lists every time wrote the same files.
- **A36. Visualiser, small things.** The rows REF and REFPAIR hold the text `TODO` and DP is always 0
  (`AssemblyVisualiser.java:337-341`). Breakends are paired in the poll order of a `PriorityQueue`
  (`:732-735`), which for breakends on the same alignment is the order of the heap; reproduced.
  `ReadViewModel.java:169-170` would exit with status 0 on an error (not reachable: `create` tests the
  same condition first). `-purple_dir` is read by nothing but the switch that keeps reads cached
  (`AssemblyConfig.java:321`); the vis files are the same with and without it. Java takes 65 s to 34
  minutes for the vis cases, the port 0.5 to 5 s.

**Crash**

- **A24. `assembly/alignment/BreakendBuilder.java:600-607`** A null `homology` is logged and then
  dereferenced (`homology.invert`) when both alignments are reverse: a NullPointerException, caught by
  the caller, which leaves the breakends formed so far. Reproduced.

**Oddity**

- **A7. `assembly/AssemblyConfig.java:341`** `AssemblyLogAllReads` is set from the flag
  `asm_build_debug`. The flag `asm_log_all_reads` (lines 148 and 423) is registered and never read.
  Reproduced.
- **A8. `assembly/types/ReadAssemblyIndices.java:32`** `isValid()` returns
  `ReadIndexStart == INVALID_INDEX`, which is inverted. No caller of it was found in the ESVEE source.
- **A15. `assembly/SequenceBuilder.java:758-790`** `findConsensusNextBases` keys a `HashMap` by `byte[]`
  (identity) and can return the first candidate seen when another reaches the required count. The
  result was analysed as independent of order, but the code is fragile. Reproduced.
- **A16. `assembly/SequenceCompare.java:211`**
  `firstIndex += checkForwards ? expectedRepeatBaseDiff : expectedRepeatBaseDiff`: both arms are the
  same. Reproduced.
- **A21. `assembly/phase/PhaseSetBuilder.java:843,846`** `candidateSupport().subList(...)` is a live
  view of the list, not a copy as the surrounding code assumes. Reproduced.
- **A25. `assembly/alignment/AlignmentFilters.java:232`** `AlignmentOrderComparator` returns -1 for
  equal keys, so it is not a consistent comparator; the resulting order then depends on the sort
  algorithm. Reproduced with a port of TimSort.
- **A28. Test utilities** pass the mapping quality as `SAMRecordSetBuilder`'s default base quality, so
  test records have MAPQ 255. The ported tests keep this.
- **A29. `assembly/ReadParseState.java:510`** `mismatchPenalty()` sums with `DoubleStream.sum()`, which
  is a compensated sum, and whose exact result also depends on the JDK (the formula changed in JDK
  17). With the SBX penalties of 0.2 and 0.3 the sum can sit exactly on a permitted-penalty threshold,
  so a plain sum, or another JDK, accepts or rejects different reads. Observed: one difference on the
  public SBX input before the port summed the same way. Reproduced.

### 4.3 Depth

**Wrong result**

- **D1. `depth/DepthTask.java:357, 382, 391`** compare `read.getMateReferenceIndex()` and
  `read.getReferenceIndex()` with `==` and `!=`. Both are `Integer`, so this is object identity, which
  holds for equal values only inside the Integer cache (-128 to 127). On a sequence with index 128 or
  above, a mate on the same sequence counts as elsewhere: no mate is awaited and ref-pair support is
  lost. Chromosomes 1-22, X and Y of the usual references are not affected (index below 25). Observed
  on synthetic data: the same reads with one contig at index 132 give `REFPAIR` summing to 89 instead
  of 7,150 over 86 variants; 79 of 86 records differ. Reproduced.
- **D2. `depth/DepthAnnotator.java:116-119`** copies each record with
  `.filters(variantContext.getFilters())`. For a record whose FILTER is "." `getFilters()` is an empty
  set, not null, so every unfiltered record is written as PASS. Observed: all 2,425 records of the
  small input. Reproduced.
- **D3. `depth/DepthAnnotator.java:362`** `x.containsPosition(variant.getStart())` tests only the
  position against the poly-G regions, not the chromosome, so a variant on any chromosome at those
  positions is dropped. Observed on a synthetic set: 42 of 577 records dropped across five chromosomes.
  Reproduced.
- **D4. `depth/DepthAnnotator.java:145-157` and `233-242`** go over `HumanChromosome` names only.
  Records on any other contig are dropped from the output without a message, while "loaded N variants"
  counts them. Observed on a synthetic set. Reproduced.
- **D5. `depth/DepthAnnotator.java:107-112`** starts a new list whenever the chromosome changes,
  replacing an earlier one: a chromosome that appears in two blocks loses the first block. Observed on
  a synthetic set. Reproduced.
- **D6. The counts depend on `-proximity_distance`.** Whether a mate is awaited depends on the slice
  end (`depth/DepthTask.java:356-357`), and the single-read path on the size of the group
  (`depth/DepthTask.java:404`). Observed on the small input: against the default, 185 of 2,425 records
  differ at 0, 15 at 500 and 92 at 100000. Reproduced.
- **D7. The cap is not a bound.** `depth/DepthTask.java:516/520` count without looking at the cap. A
  capped variant is only dropped from the list when the scan at `depth/DepthTask.java:316-338` meets
  it before finding a relevant variant; `depth/DepthTask.java:325` then lowers `MinPositionIndex` by
  one whatever the dropped variant's index. Observed at `-vaf_cap 0.05`: 223 of 4,850 sample genotypes
  end over their cap, 17 by more than half. Reproduced.

**Crash**

- **D8. `depth/UnmappedRegionDepth.java:52/62/104/108`** go over all genotypes of the record
  (`getNSamples`) but size their lists by the number of configured samples. With fewer samples
  configured than the VCF has, and `-unmap_regions`: IndexOutOfBoundsException. Observed with a
  tumor-only run on a two-sample VCF. Reproduced as a failure.
- **D9. `depth/UnmappedRegionDepth.java:113`** writes into the attribute map of a genotype that has no
  attributes; that map is htsjdk's shared unmodifiable one: UnsupportedOperationException. Observed on
  a synthetic set. Reproduced as a failure.
- **D10. `depth/DepthConfig.java:73-74` with `common/FileCommon.java:86`** Without `-sample` (or
  `-bam_file`) the tumor and reference options of the other steps are looked up, which this step does
  not register: NullPointerException. Both options are in effect required although registered as
  optional (`depth/DepthConfig.java:119`). Reproduced as a failure.
- **D11. `depth/DepthAnnotator.java:209` and `222`** With `-vcf_tag_prefix` the header lines for the
  prefixed tags are only added when the input already has REF/REFPAIR lines. On a raw assembly VCF the
  writer stops at the first record ("Key PFX_REFPAIR ... isn't defined in the VCFHeader") and leaves a
  truncated `.vcf.gz`. Reproduced as a failure; the port leaves a complete header-only file and names
  the other tag in its message.
- **D12. CRAM 3.1.** htsjdk 2.23.0 cannot read it ("CRAM version 3.1 is not supported"); it is the
  default of current samtools. Not reproduced: the port reads it.

**Oddity**

- **D13. `common/FileCommon.java:177`** Without `-output_dir` the file name starts with the text
  "null" (`null<sample>.esvee.ref_depth.vcf.gz` in the working directory). Reproduced.
- **D14. `depth/DepthConfig.java:99-103`** `-specific_chr` is registered and never read. An invalid
  `-specific_regions` logs "failed to load specific regions" and the run goes on unrestricted.
  Reproduced.
- **D15. `depth/DepthAnnotator.java:134`** The warning has one placeholder and two arguments, so it
  prints the record count where the file name was meant. Reproduced.
- **D16. `depth/UnmappedRegionDepth.java`** uses the tag REF throughout (lines 55, 107, 113, 130), not
  the tag with the `-vcf_tag_prefix`: with a prefix it reads and rewrites the input's REF values. A
  variant is added once per qualifying sample (line 68), so "setting ref depth to median for N
  variants" counts some twice. Reproduced.
- **D17. `hmf-common` `BamUtils.java:21`** describes `-bam_validation` as "STRICT (default)"; the
  default is SILENT (line 24).
- **D18. Dead code.** `depth/DepthTask.java:391-394`: the inner `if` has no effect, since both paths
  return false. `depth/DepthTask.java:156/212` (`READ_CACHE_CLEAR_COUNT`) and `229/266`
  (`readGroupTotal`, always 0 in the performance line) do nothing.
- **D19. TRACE lines.** "ref limit ... reached" (`depth/DepthTask.java:526-531`) prints intermediate
  counts that depend on `HashMap` iteration order. Not reproduced: these lines can differ; the counts
  written to the VCF do not depend on the order.

### 4.4 Caller

**Wrong result**

- **C2. `caller/VariantFilters.java:219`** For a single breakend paired with a LINE breakend,
  `breakend.lineSiteBreakend().fragmentCount(genotype)` is given the first breakend's genotype, and
  the method only reads its argument. The breakend's own fragment count is added twice instead of the
  pair's being added. Reproduced.
- **C5. `caller/CallerApplication.java:360-372`** `-manual_ref_depth` is only applied when there is a
  reference sample: a tumor-only run gets the two FORMAT header lines and no values. A reference-only
  run dies with a NullPointerException at line 368, because the tumor id is empty. Both observed.
  Reproduced, the second as a failure.
- **C6. `caller/VcfWriter.java:148` and `caller/BreakendWriter.java:108`** Output is driven by
  `HumanChromosome.values()` renamed for `-ref_genome_version`. With version 37 and a VCF whose contigs
  are named with `chr`, nothing matches and the three VCFs are written header-only, with exit code 0
  and no warning. Observed. Reproduced.
- **C11. `hmf-common` `.../sv/RepeatMaskAnnotations.java:162-166` and
  `.../region/BedFileReader.java:132-137`** A chromosome that appears again after another one starts a
  new list that replaces the first, so earlier entries are lost without a message (repeat mask, target
  regions BED). Reproduced.
- **C14. `caller/SvDataCache.java:61`** A breakend on a contig outside chr1-22/X/Y is skipped, so its
  mate stays unmatched and is dropped with only the summary warning. Reproduced.

**Nondeterminism**

- **C1. `caller/BreakendWriter.java:194`** The filters of a row are joined from a
  `HashSet<FilterType>`. Enum hash codes are identity hashes, so their order in `final_breakend.tsv`
  changes from run to run. Observed: four Java runs on the same input all differ from each other, only
  in the order of the filters within rows. Not reproduced: the port writes enum order.

**Crash**

- **C4. `caller/VcfWriter.java:223-225` with `hmf-common` `.../variant/GenotypeIds.java:89`** A sample
  id that is only a substring of the VCF's sample name is accepted by `fromVcfHeader`, but the writer
  looks the genotype up by the configured id and dies with a NullPointerException after all the work,
  leaving empty output files. Observed. Reproduced as a failure: exit 1, files left empty.
- **C12. `prep/types/DiscordantStats.java:245-255`** returns null when the file cannot be parsed, and
  `caller/VariantFilters.java:111` dereferences it. Reproduced as a failure.
- **C13. `caller/HotspotCache.java:146`** A hotspot file that cannot be read is logged and the run
  goes on without hotspots, while a line with too few fields is an uncaught
  ArrayIndexOutOfBoundsException. Reproduced, the second as a failure.
- **C16. `hmf-common` `.../genome/chromosome/HumanChromosome.java:105`** `contains("")` (or `"chr"`)
  throws NumberFormatException, because an empty string counts as numeric. Reproduced as a failure.

**Oddity**

- **C3. `caller/Variant.java:220`** `x.vcfTag().equals(x)` compares a String with the enum value and is
  never true, so `addExistingFilters` adds nothing. The filters are picked up later by
  `VariantFilters.applyExistingFilters`, except `dedup`, which that method skips on purpose.
  Reproduced.
- **C7. `caller/CallerConfig.java:141`** registers `-specific_regions`; line 82 only reads
  `-specific_chr`. The option is accepted and ignored. Observed: the run equals the standard run.
  Reproduced.
- **C8. `common/FileCommon.java:177`** Without `-output_dir` the directory is null and string
  concatenation makes file names like `null<sample>.esvee.somatic.vcf.gz` in the working directory;
  the default prep directory and input VCF get the same prefix. The same as D13. Reproduced.
- **C9. `caller/CallerApplication.java:186-187`** The germline-mode log line prints the tumor ordinal
  and id under `ref(...)` and the reference's under `tumor(...)`. Reproduced.
- **C10. `caller/CallerApplication.java:114`** "missing input files: disc-stats and frag-lengths"
  passes the two file names without placeholders, so they are not shown. The files are looked up under
  the tumor id, or under the reference id in a reference-only run, which the message does not say
  either. Reproduced.
- **C15. Unused and dead code.** `common/FilterType.java` `SGL_STRAND_BIAS` and
  `caller/FilterConstants.java` `PolyGcRegions` are never used by the caller.
  `assembly/SequenceCompare.java:476-501`: `lastRepeatSkipBases` is never set, so its branch is dead.
  `assembly/types/RepeatInfo.findDualDualRepeat` cannot find anything that
  `findMultiBaseRepeat(..., 4, ...)` has not already found.
- **C17. Test code, `esvee/src/test/.../caller/CallerTestUtils.java:40`** `new GenotypeIds(0, 10, ...)`
  with two genotypes: `VariantContext.getGenotype(10)` is null, so every test variant has no tumor
  genotype and the tests run the reference-only paths. The ported tests keep this.

## 5. Performance

One machine (macOS, Apple silicon, 20 cores), which was shared with other jobs during all runs. Time
is wall-clock; memory is peak footprint from `/usr/bin/time -l` unless RSS is named.

### Small pair, the whole chain

From the status files of the end-to-end runs (`baseline/esvee/java/small_t8/status`,
`baseline/esvee/java/small_ref/status`, `baseline/esvee/rust/small_t8/status`). The Java reference
chain ran at one thread, with `-XX:hashCode=2` from the assembly step on (section 2.5).

| Step | Java, 8 threads | Java reference, 1 thread | Rust, 8 threads |
|---|---|---|---|
| prep | 55.5 s / 43.8 GB | 241.6 s / 5.7 GB | 18.8 s / 1.9 GB |
| assembly | 17.6 s / 8.0 GB | 32.4 s / 2.7 GB | 2.8 s / 1.5 GB |
| depth | 10.4 s / 4.5 GB | 22.9 s / 2.2 GB | 1.2 s / 0.06 GB |
| caller | 11.2 s / 5.1 GB RSS | 11.8 s / 4.2 GB | 1.1 s / 0.23 GB |
| total | 94.7 s | 308.6 s | 23.9 s |

The caller has no thread option. Its time in the first column is from the step's own timing file of
that run; the Java heap limit differed between the Java columns, which is why the Java memory figures
of one step differ.

### Per step, from the step notes

**Prep** (small pair, 120 million reads in, 1.4 million records out):

| Run | Wall | User | Peak RSS |
|---|---|---|---|
| Java, 8 threads | 55-58 s | 282-288 s | 40-44 GB |
| Java, 1 thread | 234-242 s | 252-259 s | 5.8 GB |
| Rust, 8 threads | 10.2 s (22.8 s with the inputs not in the page cache) | 61 s | 2.6 GB |
| Rust, 1 thread | 66.8 s | 59 s | 1.4 GB |

On the exome pair Java took 1,100 to 1,135 s and 6.5 GB at one thread; the port 60 s and 3.9 GB at 8
threads. A profile of the port at one thread shows about 53% in htslib (inflate, record parsing, CRC)
and 13% in the low-quality check.

**Assembly** (small pair, default outputs, warm cache):

| Run | Wall | RSS | Peak footprint |
|---|---|---|---|
| Java, 8 threads | 6.5 s | 26.9 GB | 6.2 GB |
| Java, 1 thread | 30.3 s | 11.8 GB | 2.4 GB |
| Rust, 8 threads | 2.7 s | 6.6 GB | 1.47 GB |
| Rust, 1 thread | 16.8 s | 6.6 GB | 1.45 GB |

RSS is dominated in both by the BWA index image of about 5 GB, which is mapped. Other Java runs at 8
threads took up to 18 s. Stages at one thread, Rust against Java: primary assembly 5.8 against 14.1 s,
phase sets 2.2 against 6.1 s, alignment 8.6 against 9.3 s (BWA itself, the same C code).

**Depth** (small input, 2,425 variants, 4.7 million reads passed to the counting):

| Run | Threads | Wall | User | Max RSS |
|---|---|---|---|---|
| Java `-Xmx100g` | 8 | 10.43 s | 23.4 s | 4.62 GB |
| Java `-Xmx8g` | 8 | 10.66 s | 24.7 s | 1.62 GB |
| Rust, 3 runs | 8 | 0.86-1.03 s | 5.6 s | 63 MB |
| Java `-Xmx32g` | 1 | 19.58 s | 22.9 s | 3.53 GB |
| Java `-Xmx8g` | 1 | 19.75 s | 23.4 s | 1.61 GB |
| Rust, 3 runs | 1 | 5.78-5.82 s | 5.4 s | 34 MB |
| Java `-Xmx8g`, CRAM 3.0 of one chromosome | 8 | 42.9 s | 48.0 s | 2.19 GB |
| Rust, the same CRAM | 8 | 1.20 s | 7.9 s | 354 MB |

A profile of the port at one thread: about 70% inflating BGZF blocks, about 25% other htslib reading,
about 1% the counting itself. Most of the inflated data is never used, since the BAM index locates a
read only to 16 kb and a slice is about 2.4 kb on average. Java has the same cost.

**Caller** (standard run on the small input, 2,425 records):

| Run | Wall | User | Max RSS |
|---|---|---|---|
| Java `-Xmx100g` | 11.21 s | 17.0 s | 5.06 GB |
| Java `-Xmx32g`, 3 runs | 10.94-11.16 s | 17.1-17.3 s | 5.04-5.29 GB |
| Java `-Xmx4g`, 2 runs | 11.25-11.90 s | 18.6-20.2 s | 2.53-2.58 GB |
| Rust, 5 runs | 1.07-1.08 s | 3.04-3.06 s | 223-254 MB |

Nearly all of the time is loading the repeat mask (718 MB of text in a gzip file) and the two PON
files. The port went from 3.5 s to 1.07 s in four steps: reference data loaded on three threads
(3.1 s), the repeat mask parsed in blocks on four threads (2.3 s), htslib's reader (1.8 s), and
zlib-rs for inflating (1.07 s). What is left is one thread inflating the repeat mask, about 0.9 s.

### Whole genome

The pair of section 2.6. Java reference chain at one thread (`baseline/esvee/java/wgs_t1/status`,
heap limit 64 GB); Rust at 8 threads (prep from `baseline/esvee/rust/wgs_t8/status`, the later steps
of the current build from `baseline/esvee/rust/wgs_javaprep_t8/status`):

| Step | Java reference, 1 thread | Rust, 8 threads |
|---|---|---|
| prep | 8,466 s / 14.1 GB | 672 s / 5.1 GB |
| assembly | 1,679 s / 30.4 GB | 110 s / 7.2 GB |
| depth | 413 s / 12.6 GB | 60 s / 0.7 GB |
| caller | 13.5 s / 6.5 GB | 2.6 s / 0.45 GB |
| total | 10,572 s (2.9 h) | 844 s (14.1 min) |

- Java was not run at 8 threads on this input. Its output then depends on the thread count, and the
  prep step needed 44 GB at 8 threads even on the small pair.
- The Java run without `-XX:hashCode=2` took about as long (assembly 1,650 s, depth 404 s, caller
  13.8 s).
- The Rust prep step ran while the Java prep step was running, and other jobs were running too.
- The prep step appears to be limited by the disk on this input: other jobs were running, and system
  time was more than a quarter of user time.

### Whole genome, assembly memory

The first whole-genome assembly run peaked at 33.4 GB. Every one of the 29.4 million reads loaded was
still allocated at the end, about 1 kB each. Two things held them. A read and its mate refer to each
other, so counting references never freed a pair. And assemblies that were built and then rejected or
deduplicated (160,000 of 245,000) were never freed, so the reads supporting them stayed too. In Java
the garbage collector frees both.

The changes, all in `src/assembly`, measured as peak footprint at 8 threads:

| Change | Peak |
|---|---|
| before | 33.4 GB |
| when a junction group is done, free its reads that nothing else refers to, and empty the assemblies it built and discarded | 13.2 GB |
| a read whose mate is freed keeps a copy of the five mate values it reports; records keep only the attributes the step reads (MC, SA, AS, NM, UQ); a junction group frees its remaining reads once all of its assemblies are past phase set building | 9.1 GB |
| assemblies left unprocessed by phase set building drop their reads after that stage | 9.1 GB (78,000 reads alive after phase sets down to 1) |
| reads allocated from their own mimalloc heap, so the memory they held goes back to the system | 8.8 GB |
| smaller read and support-read objects (Read 512 to 368 bytes, SupportRead 280 to 200 bytes) | 7.2 GB |

| Run | Threads | Wall | Peak footprint | Max RSS |
|---|---|---|---|---|
| before | 8 | 124-149 s | 33.4 GB | 31-38 GB |
| after | 8 | 124-127 s | 7.2 GB | 12.9 GB |
| before | 1 | 883 s | 33.3 GB | 28.2 GB |
| after | 1 | 803 s | 6.7 GB | 11.8 GB |

- RSS includes the BWA index image of about 5-6 GB, which is a file mapping and not part of the
  footprint.
- Each run's `assembly.tsv`, `breakend.tsv` and raw VCF equal the first whole-genome run. All 41
  option cases of section 2.2 are still identical to Java.
- What is left at the peak, early in phase set building: 4.3 million candidate reads that assemblies
  keep for phase set building, as Java does (3.3 GB), and support reads, remote regions and assemblies
  (2.8 GB).
- `ESVEE_MEM_STATS=1` adds log lines with the footprint and live object counts at each stage. Without it
  the log is unchanged.
- The separate read heap needs the default `mimalloc` feature, which uses mimalloc 3, whose heaps can be
  used from any thread. With an older mimalloc, or with `--no-default-features`, reads come from the
  global allocator and the memory that freed reads held is returned less well.

## 6. Not ported, not tested

Not ported:
- **Depth: record validation.** htsjdk's `SAMRecord.isValid` checks behind `-bam_validation STRICT` and
  `LENIENT`. On a BAM that htsjdk finds invalid (synthetic: an unmapped read with a non-zero MAPQ) Java
  with STRICT stops with status 1 and the port carries on and writes the VCF. LENIENT only prints to
  stderr in Java and gives the same VCF. On the real BAMs STRICT changes nothing in either tool.
- **Assembly:** some of the log text of `-asm_build_debug`.

Not tested, or tested only in part:
- **Whole genome:** one Illumina tumor/normal pair (section 2.6); no whole-genome SBX or Ultima data.
- **Reference genome 37** on real data, in any step. Prep and caller were run with version 37 on genome
  38 data only (nothing matches); assembly not at all; depth on genome 38 data and a synthetic set.
- **Caller on real SBX and Ultima calls.** The caller's SBX and Ultima paths were run on Illumina data
  only. Prep and assembly were run on public SBX and Ultima reads; depth has no sequencing-type option.
- **Caller filters** `strandBias`, `lineSource`, `sbxInvZeroMotif`, `invShortFragLowVaf`, SAGA-matched
  variants, and contigs outside chr1-22/X/Y: unit tests and synthetic data only.
- **A real decoy genome image** in assembly; **sambamba** as prep's `-bamtool`; **`ReadIdTrimmer`**
  switching itself off part-way through a run.
- **SAGA** was tested with synthetic resources only, in prep and assembly.
- **Platforms.** Every run was on one macOS ARM machine. The visualiser's text positions use that
  machine's font metrics; Java on another system may place text differently.

Known differences from Java that are not faults of either:
- Header lines of an input VCF are kept as text in the depth and caller steps. htsjdk re-serialises
  them, so a header not written by htsjdk (another key order or quoting) could come out differently.
  ESVEE's own VCFs are written by htsjdk. In the depth step `GT` is passed through as text as well
  (identical for the forms tried).
- The caller does not convert genotype field `GL` to `PL` as htsjdk does; ESVEE does not write `GL`.
- In the caller, a PON end region on a contig outside chr1-22/X/Y would throw in Java
  (`ChrBaseRegion.overlaps`) and is a non-match in the port. Files ending `.bgz` are decompressed as
  gzip everywhere; Java does that for the VCF only.
- A missing SAGA index image fails at start-up in the port, and only when the first matcher is created
  in Java.

## 7. How the port was made

- **Four steps in parallel.** The four steps were ported at the same time, each by a different person
  working in a copy of a common crate skeleton. Each step was checked against the Java step's output
  for the Java step's input, so no step waited for the one before it. The four were then merged into
  this crate and the chain was run end to end (section 2.5).
- **Java first.** Before porting, each Java step was run several times at one and at 8 threads to see
  whether its output was stable. That is how the reference of section 3 was found.
- **Java collection order.** Where output depends on the iteration order of a Java collection, the
  order is emulated: `HashMap` and `HashSet` with String and Integer keys, including resize rules and
  tree bins (checked against checksums taken from a JVM over 20,000 puts and removes), and
  identity-keyed sets as they iterate under `-XX:hashCode=2`.
- **TimSort.** Java's sort is ported for the places where a comparator is inconsistent (finding A25,
  and the junction comparison for contigs that rank equally), since the result then depends on the
  algorithm.
- **Compensated summation.** `DoubleStream.sum()` is emulated where a threshold is tested on its
  result (finding A29).
- **htsjdk's writing rules.** htsjdk rewrites a record when it is changed or re-written: BAM attributes
  sorted by tag with integers in the smallest type; VCF INFO, FILTER and FORMAT keys sorted, QUAL
  re-formatted, header lines sorted and standard lines repaired. These rules are ported for the files
  each step writes.
- **The same BWA.** Java aligns through gatk-bwamem-jni 1.0.5, which wraps lh3/bwa at commit
  `cb950614` (0.7.16a-r1181) and reads the genome from an index image. The port builds the alignment
  library sources of that commit, with the non-JNI part of the wrapper, into the binary, and fills in
  the options as the Java setters do (the layout of the option structure was checked against the
  offsets the Java side uses).
- **Java's faults kept.** Where a fault in Java changes the output it is reproduced and marked with a
  comment in the code, so that a fix upstream can be followed by removing the marked code.
- **Unit tests.** The Java unit tests of each package were ported with the step: 21 for prep, 101 from
  20 classes for assembly, 2 for depth, 18 for the caller.
