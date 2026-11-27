# CHORD in Rust: record of the port

`chord-rs` ports hmftools CHORD 2.1.2: the Java part, which counts mutation contexts from PURPLE's
somatic and SV VCFs, and the embedded R part (`chord/src/main/R/chord_predict.R` with the randomForest
model `CHORD.rds`), which predicts HR deficiency. The port needs no R. Java reference:
this repository's `chord/target/chord-2.1.2-jar-with-dependencies.jar` (sources following upstream
ca04e6c39d of 2026-11-20), with R 4.5.2 and randomForest 4.7-1.2 (in a private R library). Finding ids:
CH1-CH10. Options and outputs are in [README.md](README.md).

## 1. What is ported

All three Java entry points with their options:

| Java (and R) | Rust |
|---|---|
| ChordApplication, ChordConfig, ChordConstants, ChordOutput | `main.rs`, `config.rs`, `output.rs` |
| prep: ChordDataPrep, SamplePrepTask, ChordDataWriter, MutContextCount, VcfFile, SmallVariant | `prep.rs`, `vcf_file.rs` |
| snv: SnvPrep, SnvDetails (with hmf-common SnvSigUtils) | `snv.rs` |
| indel: IndelPrep, IndelVariant, IndelDetails, IndelContext, MutationType, ContextType | `indel.rs` |
| sv: SvPrep, StructuralVariant, SvContext, SvDetails (with hmf-common VariantAltInsertCoords) | `sv.rs`, hmf-common-rs `sv` |
| hmf-common RefGenomeSource with htsjdk IndexedFastaSequenceFile | `refgenome.rs` |
| predict: ChordModel, PredictConfig | `predict/mod.rs`, `main.rs` (`chord-rs predict`) |
| chord_predict.R: feature merging, randomForest predict(type = "prob"), HR status, QC, write.table | `predict/mod.rs`, `predict/forest.rs`, `predict/rfmt.rs` |
| R's readRDS (for `-chord_model_file`) | `predict/rds.rs` |
| ChordDataPrep.main | `chord-rs prep` |

`chord-rs prep` and `chord-rs predict` stand for `java -cp chord.jar com.hartwig.hmftools.chord.prep.ChordDataPrep`
and `...predict.ChordModel`; the class names are accepted too. Nothing of CHORD is left out; it draws no
plots.

### Java and R behaviour reproduced (not faults)

- **Reading:** records in file order; PASS only (FILTER "." is not a pass) unless `-include_non_pass`;
  htsjdk's VariantContext types (SNP, INDEL, SYMBOLIC, MIXED; `isSymbolicOrSV` takes symbolic records and
  insertions of more than 150 bases); htsjdk's allele checks (a REF outside ACGTN fails); a plain or
  gzipped VCF without an index.
- **The reference:** htsjdk's errors for an unknown contig or a query past a contig's end; null (then a
  NullPointerException in CHORD) for a start before 1 or an empty range; the lengths from the `.dict`
  (CH6).
- **Contexts:** the 96 SNV bins on the C/T strand (SnvSigUtils), the ambiguous-base rule (IUPAC codes
  skipped with a warning, more than 10% of SNVs fails, compared in float), the INDEL homology and repeat
  counts (the repeat count is Java's `replaceAll` of the indel sequence as a regex, which for bases is a
  count of leading copies), the SV types (an INS is any pair one base apart, as CHORD was trained) and
  length bins.
- **Errors and partial outputs:** where Java fails the port fails at the same point, with exit status 1
  and the same files: the details files written before the failing prep, the mutation contexts file
  before a failing prediction.
- **Configuration:** hmf-common ConfigBuilder (path checks with '*' and '$' allowed, `-threads` a string
  item parsed with Integer.parseInt, unknown log4j level names rejected), sample ids from `-sample`
  (String.split) or `-sample_id_file` (ConfigUtils.loadDelimitedIdFile), the sample wildcard replaced as
  a Java replacement string (`\` and `$n` rules), output names joined as Java joins them (`<dir>//<file>`).
- **The R prediction** (section 2): the column grouping by R's regular expressions with groupColumns'
  column-count check, the fractions (`x / rowSums(x)`, NaN as 0), the forest walk, the vote fractions,
  `p_hrd` as rowSums, the QC rules and remarks, `gsub("^;|;$", ...)`, write.table's number format, the
  row names as read.delim converts them (CH3) and max.col's ties (CH1).

## 2. The prediction without R

### The forest resource

`chord-rs/resources/chord_forest.tsv` (384 KB, embedded with `include_str!`) is made once by
`scripts/export_forest.R` (in this crate) from `chord/src/main/R/CHORD.rds`:

```
Rscript chord-rs/scripts/export_forest.R chord/src/main/R/CHORD.rds chord-rs/resources/chord_forest.tsv
```

It holds what `predict.randomForest(type = "prob")` uses: the classes (BRCA1, BRCA2, none), the 12
variables in the order of `rownames(object$importance)` (C.A, C.G, T.A, T.G, del.mh.bimh.1,
del.mh.bimh.2.5, DEL_0e00_1e03_bp, DEL_1e03_1e04_bp, del.none, DUP_1e03_1e04_bp, DUP_1e04_1e05_bp,
DUP_1e05_1e06_bp), `ncat` (all numeric), and for each of the 500 trees its `ndbigtree` nodes (15,072 in
all): nodestatus, bestvar, left and right daughters, nodepred and xbestsplit. The split values are written
as C99 hexadecimal floats (`sprintf("%a")`), which hold each double exactly; the script checks that R
reads every one back to the same double. (17 significant digits would hold them too, but R's own
string-to-double conversion is not correctly rounded, so that check failed inside R for one value.) The
crate checks the file as it reads it: every interior node splits on a valid variable and points to later
nodes, every terminal node votes for a valid class.

`chord-rs` also reads `.rds` files itself (`predict/rds.rs`: XDR serialization versions 2 and 3, gzip
or uncompressed, with references, environments, ALTREP sequences), so `-chord_model_file CHORD.rds` works
without R; test `embedded_forest_equals_chord_rds` checks that the forest read from CHORD.rds equals the
embedded one.

### What the prediction reproduces

- **The forest:** each tree is walked from the root, left when `x[var] <= split` (rf.c predictClassTree);
  the probability of a class is its votes over the votes cast (`countts` normalised by its row sum), so
  every probability is k/500, computed as a double division as R does.
- **p_hrd** is rowSums of the two probabilities. R sums in long double on x86-64: for k/500 values the
  exact sum fits in long double's 64-bit significand, so it equals the double sum (which arm64 R does).
- **Numbers in the file:** write.table sets R_print.digits to 15 and writes each double with formatReal
  and EncodeReal0 (fixed or scientific, whichever is narrower). `predict/rfmt.rs` follows R's
  double-precision code (scientific(), the width rules, a three-digit exponent from 1e±100).
- **Checked against R:**
  - all 125,751 vote pairs (k1/500, k2/500, k1 + k2 <= 500) and their sum (`scripts/r_vote_pairs.R`,
    test `vote_pairs_as_r` with CHORD_R_VOTE_PAIRS): identical;
  - 1,742 doubles written by R (`scripts/r_format_fixture.R`, `testdata/r/write_table_doubles.tsv`:
    special values, every k/500, random magnitudes from 1e-320 to 1e308, random decimals): identical;
  - chord_predict.R itself against `chord-rs predict` on random synthetic samples
    (`scripts/predict_check.sh` with `scripts/random_contexts.py`): 2,000 samples (seed 1) and 20,000
    (seed 2), with the embedded forest and with CHORD.rds read directly: identical, including 5
    HR-deficient ties under the defined reference of CH1. The 20,000 cover 2,798 HR-deficient, 2,697
    HR-proficient and 14,505 undetermined samples (MSI and low-count rows on purpose).

## 3. How it was checked

**Mode:** every comparison ran the port with `-java_compat`; Java ran with `-Xmx8g` (or the pipeline's
`-Xmx12240656794`) and 1 or 2 threads as the case says (CHORD uses threads only for cohorts).

**R environment of the Java references:** a private R library, and `R_ENVIRON_USER` and `R_PROFILE_USER`
set to files of the port's work folder, so that no user `~/.Rprofile` or `~/.Renviron` is read. The profile is empty, except for the cases with ties of CH1 (`seed1.R`: `set.seed(1)`, the
defined Java reference).

**What was compared:** every file the run writes (prediction, mutation contexts, details files), byte for
byte, the set of files, and the exit status. Not compared: `CHORD.tmp.rds`, the temporary copy of the model
that Java leaves behind when R fails (CH9).

**Scripts** (kept outside this repository, with the comparison outputs): a regression runner over all
cases, the Java and Rust run scripts, generators of the renamed private copies and of the synthetic samples
(a synthetic genome and PURPLE-like VCFs whose contexts are exactly a target row, each variant's context
checked against an independent implementation of CHORD's rules), the contig-edge inputs, and the R scripts
for the formatting fixture and the vote-pair check. `scripts/export_forest.R` is in this crate.

### Inputs

- **GIAB HG008** (public): the pipeline's PURPLE VCFs of HG008-T (19,254 somatic records, 10,839 PASS:
  9,227 SNVs, 1,392 INDELs; 509 SV records, 181 SVs counted). Java's output equals the pipeline's (Hartwig
  image chord:2.1.2).
- **The private pair:** renamed copies (TUMOR) of the all-Java run's PURPLE VCFs. The port's outputs equal
  that run's CHORD outputs too.
- **The Java tests' inputs** (`chord-rs/testdata`): MINIMAL_SAMPLE, EMPTY_SAMPLE, the non-standard
  nucleotide genome.
- **Synthetic samples** (no real data): rows of the random contexts table
  made into genomes and VCFs: BRCA2-like and BRCA1-like HR deficient, MSI-like (14,420 mutations), a tie
  (TIE1), an 8-sample cohort with 5 ties and an empty sample, extras (an MNV, filtered and FILTER "."
  records, a multi-allelic SNV, single breakends, an insertion-type pair, a filtered pair, a 160-base
  insertion as sequence), plain VCFs, soft-masked genomes, no `.dict`, ids "0123" and "NA", few variants,
  SNVs only; contig-edge cases on a hand-made genome (`inputs/edge`).

### Results

`regress.sh`: **53 cases, all identical** (files and exit status). Time is wall time and memory the peak
resident set (`/usr/bin/time -l`; for Java the JVM's, R runs as a child process). Medians of 3 runs for
the first two lines, single runs for the rest.

| Case | Java exit | Java time, memory | Rust time, memory | Result |
|---|---|---|---|---|
| hg008, pipeline command | 0 | 1.17 s, 0.518 GB | 0.17 s, 0.048 GB | identical (2 files) |
| private pair, pipeline command | 0 | 0.93 s, 0.271 GB | 0.06 s, 0.022 GB | identical (2 files) |
| hg008 `-purple_dir -write_detailed_files -output_id -log_level DEBUG -threads 2` | 0 | 1.17 s, 0.499 GB | 0.18 s, 0.049 GB | identical (5 files) |
| hg008 `-include_non_pass` | 0 | 1.20 s, 0.578 GB | 0.19 s, 0.076 GB | identical |
| hg008 `-chord_model_file CHORD.rds` | 0 | 1.19 s, 0.504 GB | 0.17 s, 0.047 GB | identical |
| private pair, details, non-pass, `-log_debug` | 0 | 0.94 s, 0.281 GB | 0.07 s, 0.025 GB | identical (5 files) |
| cohort HG008-T and the private tumor, 2 threads | 0 | 1.24 s, 0.621 GB | 0.17 s, 0.059 GB | identical |
| hg008 prep only (ChordDataPrep) | 0 | 0.94 s, 0.509 GB | 0.17 s, 0.048 GB | identical (1 file) |
| hg008 SNVs and INDELs only / SVs only | 0 | 1.14 / 0.83 s | 0.17 / 0.02 s | identical |
| hg008 `-log_level WARN` | 0 | 1.17 s, 0.502 GB | 0.18 s, 0.049 GB | identical |
| hg008 `-log_level ALL`, `OFF` (CH2) | 1 | 1.39 s, 0.751 GB | 0.17 s, 0.050 GB | identical (contexts only) |
| hg008 `-log_level NOPE` | 1 | 0.34 s, 0.094 GB | 0.01 s, 0.002 GB | identical (no files) |
| Java tests: MINIMAL_SAMPLE, EMPTY_SAMPLE (details), both as a cohort, non-standard nucleotides | 0 | 0.66-0.71 s, 0.11-0.13 GB | 0.02-0.03 s, 0.005-0.007 GB | identical |
| MINIMAL_SAMPLE on a soft-masked genome (CH4) | 1 | 0.43 s, 0.108 GB | 0.01 s, 0.004 GB | identical (no files) |
| synthetic BRCA2-like, BRCA1-like (details), MSI-like | 0 | 0.68-0.81 s | 0.02-0.06 s | identical |
| synthetic tie (CH1, seed1 reference) | 0 | 0.90 s, 0.296 GB | 0.07 s, 0.025 GB | identical |
| synthetic 8-sample cohort with ties, 2 threads / 1 thread (seed1) | 0 | 1.27 / 1.24 s, 0.44-0.52 GB | 0.22 / 0.28 s, 0.10 GB | identical |
| `-sample_id_file` (3 of the cohort; seed1) | 0 | 1.14 s, 0.415 GB | 0.22 s, 0.093 GB | identical |
| `-sample C1,C1` (CH3) | 1 | 0.79 s, 0.261 GB | 0.04 s, 0.016 GB | identical (contexts only) |
| extras, extras `-include_non_pass` (details) | 0 | 0.67-0.69 s | 0.02 s | identical (5 files) |
| plain (uncompressed) VCFs | 0 | 0.68 s, 0.111 GB | 0.02 s, 0.004 GB | identical |
| soft-masked genome, INDEL-only samples (CH4) | 0 | 0.67-0.69 s | 0.02 s | identical |
| INDELs in repeats on the upper-case genome (details) | 0 | 0.69 s, 0.118 GB | 0.02 s, 0.006 GB | identical |
| no `.dict` (CH6) | 1 | 0.40 s, 0.108 GB | 0.01 s, 0.003 GB | identical (no files) |
| sample id "0123" (CH3) | 0 | 0.70 s, 0.114 GB | 0.02 s, 0.005 GB | identical ("123" in the prediction) |
| sample id "NA" (CH3) | 1 | 0.63 s, 0.116 GB | 0.01 s, 0.005 GB | identical (contexts only) |
| few variants (details), SNVs only | 0 | 0.69-0.71 s | 0.02 s | identical |
| SNV at a contig's first / last base, INDEL at its start / end (CH5) | 1 | 0.37-0.40 s | 0.01-0.02 s | identical (the INDEL-end case: the SNV details file) |
| breakend pair on an HLA contig (CH8, details) | 0 | 0.67 s, 0.113 GB | 0.02 s, 0.005 GB | identical |
| SNV on a contig missing from the reference | 1 | 0.40 s, 0.107 GB | 0.01 s, 0.003 GB | identical |
| predict only: 2,000 random samples, embedded / CHORD.rds (seed1) | 0 | 0.69 / 0.70 s, 0.11 GB | 0.06 / 0.05 s, 0.009 / 0.013 GB | identical |
| predict only, output file without a directory (CH7) | 1 | 0.33 s, 0.094 GB | 0.01 s, 0.002 GB | identical |
| configuration errors: no `-ref_genome`, `-threads x`, missing output directory, no inputs, no sample, wildcard to a missing VCF | 1 | 0.35-0.92 s | 0.01-0.17 s | identical |

**The corrected (default) mode against Java:**
- identical on HG008, the private pair and every case without a finding's trigger (the pipeline's inputs
  are not affected by any of the findings);
- different where a finding applies: the ties (synthetic tie, cohorts, sample id file: CH1), the log
  levels ALL and OFF (a prediction is written: CH2), the ids "0123", "NA" and C1,C1 (CH3), the soft-masked
  genomes (CH4), the contig edges (CH5), no `.dict` (CH6), predict without a directory (CH7), the HLA
  breakend pair (CH8). The contig missing from the reference fails in both modes.

### Unit tests

`cargo test --release`: 41 tests pass (one, `vote_pairs_as_r`, reads R's file only when CHORD_R_VOTE_PAIRS
is set). The Java tests are ported one module per class (`src/tests/`): SnvPrepTest (3), IndelPrepTest (2),
IndelVariantTest (5), SvPrepTest (3), StructuralVariantTest (3), and ChordApplicationTest (2), which Java
ignores because it needs R and which runs in the port. `src/tests/port.rs` tests each finding in both
modes (CH1, CH3, CH4, CH5, CH6, CH8), the QC rules and remarks, groupColumns' check, the embedded forest
against CHORD.rds, the HG008 prediction, the forest file checks, the cohort file, wildcards, the sample id
file and the PASS filter; `predict/rfmt.rs` tests R's number formatting (the fixture), R's id conversion,
R's Mersenne-Twister (`set.seed(1); runif(5)`) and hex floats. hmf-common-rs: 106 tests pass, among them the new
one for VariantAltInsertCoords.

## 4. Findings

Each with a switch is corrected by default behind `compat::java("CHn")`; `-java_compat` gives Java's
behaviour.

- **CH1 (chord_predict.R, max.col): the HRD type of a tie is drawn at random.**
  - **What R does:** `hrd_type <- colnames(probs_brca_type)[max.col(probs_brca_type)]` with max.col's
    default `ties.method = "random"`: probabilities within 1e-5 of the row's larger one are ties, broken
    with R's random generator, which the script never seeds (R seeds it from the time and process id).
  - **Effect:** an HR-deficient sample with p_BRCA1 = p_BRCA2 gets BRCA1_type or BRCA2_type at random.
    The probabilities are multiples of 1/500, so a tie is an exact equality: 5 of the 2,798 HR-deficient
    random synthetic samples. Rows with both probabilities 0 draw too, which shifts the draws of later
    samples in a cohort.
  - **Evidence:** synthetic sample TIE1 (p_BRCA1 = p_BRCA2 = 0.35, HR_deficient): five Java runs with
    R's default seeding gave BRCA2_type, BRCA1_type, BRCA2_type, BRCA2_type, BRCA1_type
    (five unseeded Java runs).
  - **Defined Java reference:** R with a profile calling `set.seed(1)`. In Java mode the port runs R's
    Mersenne-Twister from set.seed(1) and draws as max.col does: one uniform per tie row in row order,
    BRCA2 when 2u < 1. Identical on TIE1, the 8-sample cohort at 1 and 2 threads, the sample id file
    case and the 20,000 random samples.
  - **Corrected:** a tie goes to BRCA1_type (`ties.method = "first"`), so the output is deterministic. A
    tie has no type; "first" keeps hrd_type to its three values, which ORANGE reads.
- **CH2 (ChordModel and chord_predict.R): `-log_level ALL` or `OFF` makes the prediction fail.**
  - **What Java does:** it passes the logger's level to the R script, whose LOG_LEVEL list has TRACE to
    FATAL only. For ALL and OFF, both valid log4j levels that ConfigBuilder accepts, the first log message
    stops the script: "argument is of length zero".
  - **Effect:** exit 1 after the mutation contexts file; no prediction (and CHORD.tmp.rds is left, CH9).
  - **Evidence:** cases hg008_log_all, hg008_log_off.
  - **Corrected:** the prediction does not depend on the log level.
- **CH3 (chord_predict.R, read.delim(row.names = 1)): sample ids are converted to numbers or logicals.**
  - **What R does:** read.table converts the sample id column with type.convert before it becomes the
    row names. A column whose ids all read as logical (T, F, TRUE, FALSE), integer or double is converted,
    and the prediction file gets as.character of the values: "0123" becomes "123", "T" "TRUE", "1e5"
    "1e+05", "12.50" "12.5". "NA" stops R ("missing values in 'row.names' are not allowed"), as do two
    ids that convert to the same value or a sample given twice in `-sample` ("duplicate 'row.names'").
  - **Effect:** the prediction's sample column differs from the mutation contexts file's for digit-only
    ids with leading zeros and the like; a sample named NA, or a cohort naming a sample twice, gets no
    prediction (exit 1).
  - **Evidence:** cases synth_numeric_id (Java writes 123), synth_na_id and synth_duplicate_ids (exit 1);
    the conversions in test `row_names_as_r_converts_them` were taken from R.
  - **Corrected:** the ids as given; NA and repeated ids are predicted.
- **CH4 (RefGenomeSource / htsjdk): soft-masked reference bases are not upper-cased.**
  - **What Java does:** htsjdk returns the bases as stored. CHORD compares them with upper-case alleles
    and bin names. An SNV next to a lower-case base has an "invalid trinucleotide context" (for example
    'T>G_gtt'), which ends the run. An INDEL's lower-case flanks never match its sequence: no repeat and
    no homology, so every INDEL in soft-masked sequence is counted without context.
  - **Effect:** with a soft-masked reference (UCSC hg38 is one), CHORD fails on any sample with SNVs.
    For INDELs alone it gives silently different counts: the BRCA2-like synthetic sample (HR_deficient,
    BRCA2_type, p_hrd 0.658 on its upper-case genome) is HR_proficient with p_hrd 0 on a soft-masked copy
    of its genome. The HMF reference is upper case, so the pipeline is not affected.
  - **Evidence:** cases java_minimal_soft_masked (exit 1), synth_soft_masked (the prediction above),
    synth_soft_masked_indels (20 + 5 + 20 INDELs in repeats and 15 with microhomology all counted as
    none); test `ch4_soft_masked_reference`.
  - **Corrected:** reference bases upper-cased.
- **CH5 (SnvDetails, IndelVariant): a variant at a contig's edge ends the run.**
  - **What Java does:**
    - an SNV at a contig's first base asks for the base before it: getBases returns null, then a
      NullPointerException;
    - an SNV at the last base asks past the end: htsjdk throws;
    - an INDEL whose left flank would start before the contig (POS up to the indel length) gets null;
    - so does one whose right flank would start after its end (an INDEL at the last bases).
  - **Effect:** exit 1. Human contigs start and end in N runs where nothing is called. Circular chrM, or
    small contigs, could have such variants; HG008 and the private pair have none.
  - **Evidence:** cases edge_snv_first, edge_snv_last, edge_indel_start and edge_indel_end (Java writes
    the SNV details file, then stops at the INDEL); test `ch5_contig_edges`.
  - **Corrected:** for SNVs, positions outside the contig are read as 'N'. The context is then ambiguous,
    so the SNV is skipped with a warning and counts towards the 10% limit. INDEL flanks are cut at the
    contig ends.
- **CH6 (RefGenomeSource.getChromosomeLength): chromosome lengths need a `.dict` file.**
  - **What Java does:** the length comes from htsjdk's sequence dictionary (`<name>.dict` or
    `<file>.dict`), which htsjdk does not build from the `.fai`. Without one, the first INDEL ends the run
    with a NullPointerException.
  - **Evidence:** case synth_no_dict; the Java tests' dummy genome without its `.dict`; test
    `ch6_no_dictionary`.
  - **Corrected:** lengths from the `.fai`. In Java mode the `.dict` lengths are used.
- **CH7 (ChordModel.main): predicting into a file without a directory fails.**
  - **What Java does:** the embedded model is first copied to `new File(output).getParent() +
    "/CHORD.tmp.rds"`. getParent() is null for a bare file name, so the copy goes to
    "null/CHORD.tmp.rds" and fails, and the run exits 1.
  - **Evidence:** case predict_no_dir.
  - **Corrected:** no temporary file; the predictions are written.
- **CH8 (VariantAltInsertCoords.fromRefAlt): a breakend's contig name is cut at its first ':'.**
  - **What Java does:** the other breakend's "chromosome:position" is split at every ':' and its first
    two parts taken. "HLA-A*01:01:01:01:2100" becomes chromosome "HLA-A*01", position 1, so a breakend
    pair within one such contig is counted as a translocation.
  - **Effect:** contigs with ':' in their names only (the HLA contigs of GRCh38_masked_exclusions_alts_hlas);
    none in HG008 or the private pair.
  - **Evidence:** case edge_hla_sv: Java counts TRA (details: BND, other breakend "HLA-A*01" at 1),
    corrected DEL of 2,000 bases; test `ch8_contig_names_with_colons`.
  - **Corrected:** split at the last ':'.
- **CH9 (ChordModel): the temporary copy of the model.** No switch.
  - **What Java does:** each run copies CHORD.rds to `<output_dir>/CHORD.tmp.rds` and deletes it in a
    finally block. When R fails, the catch block calls System.exit, which skips the finally block, so the
    119 KB file is left (cases hg008_log_all, hg008_log_off, synth_na_id, synth_duplicate_ids). Two runs
    sharing an output directory share the temporary path: each overwrites it and the first to finish
    deletes it under the other.
  - **The port** has no temporary file. The comparisons leave Java's leftover out.
- **CH10 (chord_predict.R logMessage): the R log's time shows the hour twice.** No switch; log only.
  `format(Sys.time(), "%H:%H:%OS3")`: a run at 03:29 logs 03:03:31.093. The port writes its log lines in
  the Java layout.

Kept as Java does, as the intent is unclear or it is the model's definition:
- an INS is any breakend pair one base apart, whatever the inserted bases (CHORD was trained so, says
  the code); INS and single breakends are not counted;
- complex INDELs (both alleles longer than one base): the "indel sequence" is the longer allele less its
  first base, longer than the length difference that sets the flanks; SAGE emits few;
- breakends are taken once per pair by ID and MATEID, so records with ID "." after the first are skipped;
- MNVs are not counted; a multi-allelic record counts its first ALT;
- in a cohort, TaskExecutor's result is ignored: a failed sample ends the run when the file is written
  (NullPointerException), and for a small cohort the file is left empty, as the buffered writer is never
  flushed. The port exits 1 with an empty file (it does not reproduce larger partial files);
- SnvPrepTest and IndelPrepTest's `providingWrongVcfTypeThrowsError` expect an IllegalStateException that
  comes from the missing reference genome, not from the VCF type; ported as they are;
- SvPrep logs "Writing indel details to" for the SV details file (log only).

## 5. Speed and memory

At the pipeline's command the port is 7 times faster on HG008 (0.17 s against 1.17 s) and 15 times on the
private pair (0.06 s against 0.93 s), in a tenth of the memory (0.048 against 0.52 GB; 0.022 against
0.27 GB). Predicting 20,000 samples takes the port 0.39 s and 61 MB, R 0.99 s and 205 MB.

Where the time goes (HG008, `-log_level DEBUG`, 0.16 s in the process):
- the somatic VCF (3.6 MB bgzipped, 19,254 records): 0.10 s, inflating and parsing, once (Java reads it
  twice, for SNVs and for INDELs);
- SNV contexts: 0.02 s, opening the reference (index and dictionary of 3,366 contigs) and 9,227 reads;
- INDELs 5 ms, the SV VCF 3 ms, the prediction 4 ms (2 ms of it reading the embedded forest);
- the remaining 0.02 s pass between the end of SvPrep and the writing of the counts, where the records
  and the reference are dropped.

Java's 1.2 s are the JVM, htsjdk reading the VCF twice, and R: starting Rscript, loading randomForest,
readRDS of the model and the prediction take about 0.3 s.

The SV VCF first took 0.33 s to read: hmf-common-rs's `VcfHeader::from_lines` checks each `##contig` line
against all earlier lines, which is quadratic in the 3,366 contig lines of PURPLE's GRCh38 headers. chord-rs
reads the header itself and builds the VcfHeader from the INFO, FORMAT and FILTER lines only, which is all
the record parsing needs (3 ms). hmf-common-rs is unchanged; a set of contig ids there would save the same
0.3 s per VCF in the other ports (PURPLE, LINX).

**Threads:** CHORD uses threads only for cohorts, one sample per thread; the 8-sample cohort takes 0.28 s
at 1 thread and 0.22 s at 2.

## 6. Open points

- **GRCh37** is not tested (no 37 inputs at hand); CHORD has no version-specific code beyond the
  reference.
- **CH3 in Java mode** reproduces R's conversion of ids that read as logicals, integers or decimal
  numbers; hexadecimal and complex numbers, and integer-valued doubles of more than 15 digits (which R
  4.3 and later print in full), are not converted as R would.
- **R on x86-64** formats numbers through long double; the port follows R's double code (arm64), which
  was the R checked here. For vote fractions the two agree (each value is within a few ulps of a short
  decimal); for a custom model of more than 1,024 trees, p_hrd could differ in the last bit on x86-64 R.
- **`-chord_model_file`:** RDS files compressed with bzip2 or xz are rejected (save with
  `compress = "gzip"`); forests with categorical predictors are not supported (CHORD's are numeric).
- **Malformed VCFs:** alleles are checked as htsjdk does; other htsjdk record checks are not all
  reproduced, and error texts differ.
- **oncoanalyser:** the CHORD module's versions.yml runs `R --version` and `Rscript -e
  'packageVersion("randomForest")'`, which fail without R; running the port needs that part of the
  module changed, as for the other ports (`--rust_ports_bin`).

## 7. Shared code

Added to `hmf-common-rs` (nothing changed):
- **`sv::VariantAltInsertCoords`** with `from_ref_alt` (hmf-common VariantAltInsertCoords.fromRefAlt:
  orientation, inserted sequence and other breakend of an SV ALT, with Java's substring errors) and
  `from_ref_alt_split` (the other breakend's position after the last ':', for CH8), and a test.

Used from hmf-common-rs: `config_builder`, `compat`, `log`, `refgenome`, `variant_context` (records and
types), `sv` (StructuralVariantType, MATE_ID, match_breakend_alt), `gz` (line reading, gzip for RDS),
`java_text`, `fmt`.

## In the pipeline (GIAB HG008, 2026-11-27)

HG008 rerun 9 on Google Batch ran SIGS, Virus Interpreter, PEACH and CHORD as ports, in Java mode, with every
other step from the first run's cache. CHORD's outputs (both files) are identical to the first run's Java
outputs. CUPPA and ORANGE, which ran again on them, gave the same outputs as before. CHORD took 7.4 s, 51 MB (mostly staging: 6% CPU),
against Java's 9 s, 439 MB. The rerun's tasks cost $0.02.
