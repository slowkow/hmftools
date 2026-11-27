# Virus Interpreter in Rust: record of the port

This report covers how `virus-interpreter-rs` was checked against Java Virus Interpreter 1.7.1, the faults found in the Java code, and the
timings. Java Virus Interpreter 1.7.1 is this repository's `virus-interpreter/target/virus-interpreter-1.7.1-jar-with-dependencies.jar`, built from sources that follow upstream commit ca04e6c39d
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

HG008's BamTools metrics were first replaced by a synthetic summary (mean coverage 61.25), as only the
pipeline's bucket had them. HG008's VIRUSBreakend summary is empty (no viruses), so that case did not
depend on them. A case with the real summary (`hg008_real`) was added at the merge.

**Against the pipeline's own outputs** (`virus-interpreter:1.7.2-beta.1`): for HG008 and the private pair, the jar's outputs equal the
pipeline's, file for file.

## What is ported

All of `VirusInterpreterApplication`.

| Java | Rust |
|---|---|
| VirusInterpreterApplication, VirusInterpreterConfig | `main.rs` |
| VirusInterpreterAlgo, CoveragesAnalyzer, CoveragesAnalysis | `algo.rs` |
| TaxonomyDb(File), VirusReportingDb(File, Model), VirusBlacklistingDbFile; hmf-common VirusBreakendFile, VirusBreakend, VirusBreakendQCStatus, AnnotatedVirusFile, AnnotatedVirus, VirusType, VirusLikelihoodType | `files.rs` |
| hmf-common PurpleQCFile (the QC statuses) | hmf-common-rs `purple_files::read_purple_qc` (added) |
| hmf-common PurityContextFile (purity, ploidy), BamMetricSummary | hmf-common-rs `purple_files`, `metrics` (existing) |

## The option name of the BamTools folder

The jar (1.7.1, the fork's sources) registers `-bam_metrics_tumor_dir`. Oncoanalyser passes
`-tumor_metrics_dir`, which Hartwig's 1.7.2-beta.1 image takes; the jar rejects it ("unregistered config
item"). The port accepts both, so it runs with the pipeline's command line. The Java runs used the jar's name,
and the port runs of the pipeline cases the pipeline's.

## Java behaviour reproduced

- **The databases:** each kept as Java keeps it:
  - the taxonomy database is read without a header line, and a line counts only with exactly two
    tab-separated fields once trailing empty ones are dropped;
  - the reporting database needs six fields, the blacklist four; other lines are warned about;
  - a later line of a taxid replaces an earlier one;
  - fields are trimmed as `String.trim` trims them, but the driver likelihood is not trimmed;
  - a taxid that is not a number, a virus name or likelihood that is not known, or text that is not UTF-8
    ends the run.
- **The summary:** VIRUSBreakend's summary is split with trailing empty fields dropped, so an empty QC status
  (the last column) reads as NO_ABNORMALITIES; a line of fewer than 24 fields ends the run.
- **The output:** doubles as Java's `Double.toString` writes them; `null` for no interpretation and for the
  expected coverage under failing PURPLE QC; the likelihood by its enum name.

## Cases and results

| Case | Java time, memory | Rust time, memory | Exit | Result |
|---|---|---|---|---|
| hg008 (pipeline command; VIRUSBreakend found no virus) | 0.82 s, 938 MB | 0.13 s, 6.1 MB | 0/0 | identical |
| t0 (pipeline command; no virus) | 0.79 s, 899 MB | 0.13 s, 6.1 MB | 0/0 | identical |
| 13 synthetic viruses, private pair's PURPLE and metrics | 0.79 s, 923 MB | 0.16 s, 6.1 MB | 0/0 | identical |
| the same on HG008's PURPLE files | 0.82 s, 911 MB | 0.16 s, 6.0 MB | 0/0 | identical |
| the same with PURPLE QC FAIL_NO_TUMOR, `-log_debug` | 0.83 s, 923 MB | 0.15 s, 6.1 MB | 0/0 | identical |
| the same with FAIL_CONTAMINATION | 0.82 s, 925 MB | 0.16 s, 6.1 MB | 0/0 | identical |
| a summary with a header only | 0.80 s, 945 MB | 0.13 s, 6.0 MB | 0/0 | identical |
| a reported virus blacklisted (VI2) | 0.79 s, 923 MB | 0.15 s, 6.1 MB | 1/1 | same failure |
| no `-output_dir` (VI1) | 0.80 s, 935 MB | 0.15 s, 6.1 MB | 0/0 | identical |
| no metrics file | 0.82 s, 954 MB | 0.15 s, 6.1 MB | 1/1 | same failure |

The **13 synthetic viruses** cover:
- EBV non-integrated above and below the expected coverage, and integrated below 90%;
- EBV's parent taxid;
- HPV 16 integrated (no minimum) and not integrated (minimum not met);
- a LOW-likelihood HPV, an HPV failing VIRUSBreakend's QC (LOW_VIRAL_COVERAGE), HBV;
- a species outside the reporting database;
- a blacklisted genus (HIV);
- a reference taxid missing from the taxonomy (an empty name);
- HHV-8 with no coverage;
- every QC status of VIRUSBreakend.

**Corrected mode:** identical to Java in all cases but no `-output_dir` (VI1: the file's name) and the error
message of VI2.

## Findings

- **VI1 (parseOutputDir with no `-output_dir`):** as SI5; Java writes `null<sample>.virus.annotated.tsv` in the
  working directory. Corrected: the working directory.
- **VI2 (VirusInterpreterAlgo.analyze): the message of the "reported and blacklisted" error.** Java writes
  "Virus with taxid {} configured as reported and as blacklisted10376": the placeholder is never filled
  and the taxid is appended. Corrected: "virus with taxid 10376 configured as reported and as blacklisted".
  Both exit 1. This happens only if the reporting and blacklisting databases overlap (they do not).
- **VI3 (TaxonomyDbFile.loadFromTsv): the whole taxonomy database is held in memory.** No switch; memory only.
  - **What Java does:** it keeps all 2.58 million names of the 89 MB file in a HashMap, and peaks at about 0.9 GB.
  - **What the port does:** it reads the summary first and keeps the names of its reference taxids only. It
    still parses every line's taxid and checks its fields and UTF-8, so a file Java rejects is rejected.

Kept as Java does (intent unclear): PURPLE quality is unacceptable only for FAIL_NO_TUMOR and
FAIL_CONTAMINATION. The newer FAIL_TINC does not count.

## Speed and memory

Five to six times faster, in 6 MB against about 920 MB. Nearly all of the 0.13 s is the scan of the taxonomy
database (89 MB, 2.58 million lines, read in 4 MB blocks and each line's taxid parsed). The other files are a
few lines each.

## Open points

- No case with a real virus: both whole-genome pairs have empty VIRUSBreakend summaries, so the virus cases
  are synthetic (real taxids of the resource databases).
- HG008's BamTools summary is not local (see section 0).

## Shared code

Added to `hmf-common-rs` for this port (additions only; no existing behaviour changed):

- `purple_files`: `PURPLE_QC_STATUSES`, `PurpleQc` (with `has_status`), `purple_qc_filename` and `read_purple_qc`.
  This is PurpleQCFile.read with Java's defaults, its old FAIL_ names read as WARN_ statuses, and its errors
  for unknown names or a missing gender line. Test `purple_qc_tests::reads_the_qc_file`. Virus Interpreter
  reads the QC statuses with it; ORANGE will need it too.
