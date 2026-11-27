# PEACH in Rust: record of the port

This report covers how `peach-rs` was checked against Java PEACH 2.0, the faults found in the Java code, and the
timings. Java PEACH 2.0 is this repository's `peach/target/peach-2.0-jar-with-dependencies.jar`, built from sources that follow upstream commit ca04e6c39d
of 2026-11-20. Options and outputs are in [README.md](README.md).

## How it was checked

The port ran with `-java_compat` against the jar (Java 23, `-Xmx8g`; the program has no thread option) on the
same inputs, one run at a time on a shared Mac. "Identical" means every file Java wrote equals the port's byte
for byte, the same files exist, and the exit codes agree. Where Java fails, the port fails on the same input
and leaves the same files written. Times are wall times (a few tens of milliseconds are noise at a load
average of about 4), and memory is the peak resident set (`/usr/bin/time -l`).

Inputs:
- **GIAB HG008** (public), from the oncoanalyser run of 2026-11-26, with the pipeline's command line and
  resources.
- **A private whole-genome pair** from an all-Java oncoanalyser run, as renamed copies (TUMOR / NORMAL).
- **Synthetic inputs** made for the port, with no sample data in them.

**Against the pipeline's own outputs** (PEACH from its module's container): for HG008 and the private pair, the jar's outputs equal the
pipeline's, file for file.

## What is ported

All of `PeachApplication`.

| Java | Rust |
|---|---|
| PeachApplication, PeachConfig, PeachUtils, PeachQCStatus | `main.rs`, `algo.rs` |
| PeachAlgo, HaplotypeCaller, HaplotypeAnalysis, HaplotypeCombination | `algo.rs` |
| event (HaplotypeEvent, VariantHaplotypeEvent, HaplotypeEventFactory), haplotype (DefaultHaplotype, NonDefaultHaplotype), panel (GeneHaplotypePanel, HaplotypePanel), PanelLoader; HumanChromosome.fromString | `panel.rs` |
| HaplotypeEventLoader, DrugInfoLoader, HaplotypeFunctionLoader, DrugInfoStore, HaplotypeFunctionStore | `loaders.rs` |
| output (EventsFile, EventsPerGeneFile, AllHaplotypeCombinationsFile, PeachGenotypeExtractor, QcStatusFile); hmf-common PeachGenotypeFile | `output.rs` |

## Java behaviour reproduced

- **No order dependence:** every output is sorted, by `String.compareTo` and the combination comparator, so
  nothing depends on Java's HashMap order; the port keeps its maps in BTreeMaps.
- **Errors where Java fails:** a duplicated event, a relevant multi-allelic record, more than one function or
  URL for a haplotype or drug of a called gene, a sample not in the VCF, and panel errors (one default and
  one wild type per gene, events without a matching id) all end the run. A failure while writing leaves the
  files Java leaves; the duplicated-drug case writes three of the five.
- **Event ids:** `VAR_chr<n>_<pos>_<ref>_<alt>` whatever the VCF's chromosome names.

## Cases and results

| Case | Java time, memory | Rust time, memory | Exit | Result |
|---|---|---|---|---|
| hg008 (pipeline command, the normal) | 0.39 s, 143 MB | 0.02 s, 2.3 MB | 0/0 | identical (5 files) |
| t0 (pipeline command) | 0.40 s, 144 MB | 0.01 s, 2.3 MB | 0/0 | identical (5) |
| hg008, the tumor sample | 0.40 s, 143 MB | 0.02 s, 2.3 MB | 0/0 | identical (5) |
| t0, the tumor sample | 0.37 s, 141 MB | 0.02 s, 2.3 MB | 0/0 | identical (5) |
| hg008, panel only (no drugs or functions) | 0.40 s, 145 MB | 0.01 s, 2.3 MB | 0/0 | identical (5) |
| the Java tests' VCF, panel, drugs and functions | 0.30 s, 87 MB | 0.01 s, 2.2 MB | 0/0 | identical (5) |
| the Java tests' plain VCF, overlapping-genes panel | 0.30 s, 85 MB | 0.00 s, 1.9 MB | 0/0 | identical (5) |
| synthetic: DPYD *1/*2A, UGT1A1 *28/*28 | 0.28 s, 89 MB | 0.00 s, 2.1 MB | 0/0 | identical (5) |
| the same as a plain VCF | 0.28 s, 86 MB | 0.00 s, 2.0 MB | 0/0 | identical (5) |
| half of *B3; a no-call at *28 (FAIL_NO_COMBINATION_FOUND, FAIL_EVENT_WITH_UNKNOWN_COUNT) | 0.29 s, 89 MB | 0.00 s, 2.2 MB | 0/0 | identical (5) |
| *2A and phased *B3; *28 and *6 homozygous (WARN_TOO_MANY_ALLELES_FOUND) | 0.29 s, 90 MB | 0.01 s, 2.2 MB | 0/0 | identical (5) |
| a filtered *2A, another allele at its position, a deletion over *28's site | 0.30 s, 90 MB | 0.00 s, 2.1 MB | 0/0 | identical (5) |
| 2846A>T/*13 with hom-ref events, *36/*37 | 0.29 s, 89 MB | 0.00 s, 2.1 MB | 0/0 | identical (5) |
| header only | 0.27 s, 86 MB | 0.01 s, 2.1 MB | 0/0 | identical (5) |
| an output folder that does not exist | 0.29 s, 90 MB | 0.01 s, 2.1 MB | 0/0 | identical (5) |
| an irrelevant multi-allelic record (PE1) | 0.27 s, 84 MB | 0.01 s, 2.1 MB | 1/1 | same failure |
| a record on chrM (PE1) | 0.27 s, 86 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| an irrelevant mixed genotype `./1` (PE1) | 0.29 s, 86 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| an irrelevant record without ALT (PE1) | 0.27 s, 86 MB | 0.01 s, 2.1 MB | 1/1 | same failure |
| a multi-allelic record at a panel position | 0.27 s, 86 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| an event twice | 0.27 s, 86 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| a sample not in the VCF | 0.28 s, 88 MB | 0.00 s, 2.1 MB | 1/1 | same failure |
| a drug listed twice | 0.30 s, 90 MB | 0.00 s, 2.2 MB | 1/1 | same failure, 3 files identical |
| a function file without its `function` column (PE2) | 0.25 s, 80 MB | 0.00 s, 1.9 MB | 1/1 | same failure |
| a function listed twice for a haplotype not called | 0.29 s, 90 MB | 0.00 s, 2.1 MB | 0/0 | identical (5) |

**Corrected mode:** identical to Java in all but the four PE1 cases, where the port completes. In each it
calls DPYD *1/*2A and UGT1A1 *1/*1, which is what the relevant records give. PE2 changes only an error message.

## Findings

- **PE1 (HaplotypeEventLoader.handleVariantContext): every record of the VCF must be a valid event, relevant
  or not.**
  - **What Java does:** it turns each record into an event and reads the sample's genotype type *before*
    checking whether the record covers a panel position. So any record anywhere in the VCF ends the run if
    it is:
    - on another contig than 1-22, X and Y (HumanChromosome.valueOf);
    - multi-allelic;
    - without an ALT;
    - or with a mixed or unavailable genotype.
  - **Evidence:** the four PE1 cases; each record is far from DPYD and UGT1A1.
  - **Count:** PURPLE's germline VCFs of HG008 and the private pair have none of these, so no effect there.
    Another caller's VCF, or one with chrM or alt contigs, would fail.
  - **Corrected:** records on other contigs or that cover no panel position are skipped first; relevant
    records keep Java's checks (a relevant multi-allelic record still ends the run). A sample missing from
    the VCF is reported by name up front (Java: a NullPointerException on the first record).
- **PE2 (HaplotypeFunctionLoader.loadFunctions): the error for a function file names the drugs file.**
  Java says "failed to load drugs TSV: <function file>". Corrected: "failed to load haplotype functions TSV".
  Message only.

Kept as Java does:
- a hemizygous call (`1`) counts as two copies (HOM_VAR); the panel genes are autosomal;
- an event at a panel position with an allele no haplotype has makes the gene unresolvable
  (FAIL_NO_COMBINATION_FOUND), as designed.

## Speed and memory

The runs take 10-20 ms of work (the log's span; wall times of 0.00-0.02 s), in about 2 MB, against Java's
0.3-0.4 s and 85-145 MB. The time is reading the germline VCF (846 records for HG008). In the corrected mode
only the records at panel positions are parsed.

## Open points

- GRCh37 data was not run beyond the Java tests' 37-coordinate panel.

## Shared code

Added to `hmf-common-rs` for this port (additions only; no existing behaviour changed):

- `variant_context::VcfReader::next_line()`: the next data line as read, unparsed, for readers that check a
  few columns before `VcfRecord::parse` (SIGS skips filtered records this way; PEACH skips records off the
  panel).

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 9 on Google Batch ran SIGS, Virus Interpreter, PEACH and CHORD as ports, in Java mode, with every
other step from the first run's cache. PEACH's outputs (all 5 files) are identical to the first run's Java
outputs. CUPPA and ORANGE, which ran again on them, gave the same outputs as before. PEACH took 616 ms, 2.8 MB,
against Java's 3 s, 80 MB. The rerun's tasks cost $0.02.
