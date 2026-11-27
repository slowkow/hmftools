# NEO in Rust: record of the port

This report covers how `neo-rs` was checked against Java NEO 1.4, the faults found in the Java code, and the timings.
Java NEO 1.4 is this repository's `neo/target/neo-1.4-jar-with-dependencies.jar`, built from sources that follow
upstream commit ca04e6c39d of 2026-11-20. Finding ids: NE. Options, outputs and tests are in [README.md](README.md).

oncoanalyser runs two steps: NEO_FINDER (`com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder`, the jar's main class)
and NEO_SCORER (`com.hartwig.hmftools.neo.score.NeoScorer`). Both are ported. RNA inputs (ISOFOX) are optional and
absent in the DNA-only runs.

## How it was checked

- **Java:** the jar, local JDK, `-Xmx8g -XX:+UnlockExperimentalVMOptions -XX:hashCode=2`, one thread (the pipeline
  passes no `-threads`), one Java run at a time.
  - NEO's HashSet orders that reach the output are of strings and integers (the alleles, a peptide's neo-epitope ids,
    the shared upstream groups), so they do not depend on identity hashes. The hashCode flag is kept for uniformity.
  - Java's output depends on thread timing only with several samples and `-threads` > 1 (NE10). The defined
    reference is `-threads 1`.
- **The port:** `-java_compat`, timed with `/usr/bin/time -l`.
- **Comparison:** every output file, byte for byte, plus the exit codes.
- **Scripts:** a run script per tool and mode (the default mode, and `-java_compat_except` for each finding), the
  input generator for the small cases and a case runner; kept outside this repository, with the outputs.

### Inputs

- **HG008** (public): the pipeline's exact NEO_FINDER and NEO_SCORER command lines (`pipeline_full/commands/`), with
  the staged names mapped to `data/resources` (`misc/neo/binding`, `misc/neo/tpm_cohort/...transcript_medians.csv.gz`,
  `common/ensembl_data`). The scorer reads the Java finder's output.
- **The private pair:** renamed copies (TUMOR / NORMAL) of:
  - PURPLE's purity, qc and somatic VCF; the VCF's sample columns are renamed too (no other line carries an id);
  - LINX's neoepitope file and LILAC's allele file;
  - the pipeline's NEO outputs.

  The same command lines as HG008. No real id is in any copied file or output (checked with grep for both ids).
- **Small cases** (HG008-based):
  - finder:
    - no fusions; no point mutations; chr1-5 only;
    - `-gene_id_file` with every fusion gene, and without one (NE3);
    - `-req_amino_acids` 10 and 20; `-write_trans_data`;
    - a sample id file with two samples at `-threads 2`;
    - a synthetic missense in an IG V-segment gene (NE14);
  - scorer:
    - no fusions; no point mutations; no neo-epitopes;
    - `-cancer_type` Lung and an unknown type;
    - `-write_types` (with the types that write nothing) and ALLELE_PEPTIDE alone;
    - thresholds and `-output_id`;
    - every score file given by option, without `-score_file_dir`;
    - only the position weights and score distribution;
    - `-score_file_dir` with a file option (NE8);
    - score files without the optional flank weights and expression likelihood distribution;
    - score files without one allele's score distribution (NE13);
    - `-check_self_recog`;
    - no TPM medians; no Ensembl data;
    - LILAC with a lost and a disrupted allele; LILAC without alleles;
    - `-sample_data_dir` with and without `-rna_somatic_vcf` (NE11);
    - a sample id file with two samples at 1 and 2 threads (NE10);
    - no `-sample`; missing PURPLE; missing neo-epitope file;
    - `-peptide_lengths` (NE9);
    - a cohort expression matrix; `-threads 2`;
  - scorer with **synthetic RNA** (`make_small_inputs.sh`; values computed from transcript names and positions):
    - ISOFOX transcript data and neo-epitope coverage, and an RNA genotype (DP, RC_CNT) added to the HG008 somatic
      VCF: with all transcripts; with a third of them missing; with a cancer type;
    - ISOFOX only; the RNA VCF only; RNA plus the cohort matrix;
    - RNA without the ISOFOX files;
    - RNA from a sample id file (NE12).

### Results (`-java_compat`)

| Case | Files | Java time, peak RSS | Rust time, peak RSS | Result |
|---|---|---|---|---|
| HG008 finder, pipeline command | 1 | 2.01 s, 1.21 GB | 0.72 s, 153 MB | identical |
| HG008 scorer, pipeline command | 2 | 8.51 s, 1.43 GB | 0.78 s, 175 MB | identical |
| Private pair finder | 1 | 1.65 s, 1.27 GB | 0.67 s, 153 MB | identical |
| Private pair scorer | 2 | 5.15 s, 1.44 GB | 0.56 s, 145 MB | identical |
| 10 small finder cases | 1-2 each | 1.41-1.76 s, 0.62-1.38 GB | 0.59-0.72 s, 7-153 MB | all identical; one exits 1 in both (NE3) |
| 36 small scorer cases | 2 each | 1.31-13.09 s, 0.54-1.66 GB | 0.20-1.16 s, 114-207 MB | 35 identical, 5 exit 1 in both; `sc_cohort_threads2` differs (NE10) |

- **Java against the pipeline:** the jar's outputs equal the pipeline's (Hartwig's neo:1.4-beta.3 image) for both
  steps, on HG008 and on the private pair.
- **Exit 1 in both:** `fin_gene_ids_nofus` (NE3), `sc_dir_and_file` (NE8), `sc_sample_dir_only` (NE11),
  `sc_missing_purple`, `sc_missing_neo`, `sc_rna_no_isofox`. Java leaves its output files empty; so does the port
  (it writes a file when the sample is done).
- **`sc_cohort_threads2`:** Java at 2 threads writes SAMPLE2's rows first. Sorted, its files equal those at 1 thread,
  and the port's files (at 2 threads) equal Java's at 1 thread.
- **Final regression:** all cases rerun with the final binary (port only): identical as above.

### The default mode

Each case was also run with `-java_compat_except NEn` for each finding, which shows which findings change it
Every difference between the default mode and Java comes from these:

| Finding | Changes |
|---|---|
| NE1, NE2, NE5 | every finder case with point mutations. HG008: NE1 6 rows (CodingBaseCigarUp), NE2 4 rows (CigarUp), NE5 3 rows (CodingBasesUp), 12 rows in all (session 1's "10 rows, CigarUp only" missed NE5). Private pair: 1 CigarUp and 1 CodingBasesUp row. |
| NE3 | `fin_gene_ids_nofus` runs instead of exiting 1 |
| NE6 | every scorer case with TPM medians: HG008 all 98 neo-epitope rows and all 64,344 peptide rows; private pair 44 of 44 and 32,280 of 32,280 |
| NE8 | `sc_dir_and_file` runs (as without the file option) |
| NE9 | `sc_peptide_lengths` scores lengths 9-10 only |
| NE11 | `sc_sample_dir_only` runs |
| NE13 | `sc_no_rank_dist`: the 10,724 A0201 rows |
| NE14 | `fin_ig_v`: the IG V-segment neo-epitope is dropped |

NE4, NE7, NE10 and NE12 change none of the cases in the default mode (NE10 and NE12 keep Java's behaviour there too).

## Findings

Each is corrected by default and reproduced by `-java_compat` (or kept as Java, where said).

- **NE1 (finder).** EpitopeUtils' upstream excerpt has `lastMatchBase - baseEnd -  - 1` (a stray minus, so `+ 1`).
  - Effect: each intron in a forward-orientation upstream cigar is 2 bases too long (`epitope.rs`).
  - HG008: 6 rows; private pair: 1.
- **NE2 (finder).** For a point mutation on an exon's first or last base, Java joins the lower excerpt, the mutated
  bases and the upper excerpt without the intron between them.
  - Example: a cigar `93M` where `45M2306N48M` is right. HG008: 4 rows.
  - Corrected: an N element for each gap.
- **NE3 (finder).** A LINX fusion whose gene is outside `-gene_id_file` (so not in the cache): Java's null
  transcript list throws a NullPointerException (exit 1).
  - Shown by `fin_gene_ids_nofus`. LINX's candidates include fusions that give no neo-epitope, so a list made
    from NEO's output alone crashes too.
  - Corrected: the fusion is skipped.
- **NE4 (finder).** A non-coding upstream transcript in a fusion: Java unboxes its null coding start
  (NullPointerException). Corrected: skipped. No case reaches it.
- **NE5 (finder).** Frameshift coding bases: Java keeps the bases of the first, unrestricted fetch (every base to the
  transcript's end), while the positions and cigar come from the restricted one.
  - On the reverse strand the extra bases come first, so a reader of the file (ISOFOX) aligns them to the wrong
    positions.
  - Corrected: the restricted fetch's bases, as SvNeoEpitope takes them. HG008: 3 rows; private pair: 1.
- **NE6 (scorer).** NeoEpitopeData.setExpressionData starts the cohort TPM sums (cancer type and pan-cancer, per
  stream) at NO_TPM_VALUE = -1 and adds each transcript's median to it, so every cohort TPM is 1 too low.
  - In a DNA-only run the cohort TPM is the TPM, so the error flows into ExpectedTpm, EffectiveTpm, the expression
    likelihood and its rank. Negative TPMs reach the expression likelihood as the lowest bucket.
  - HG008: 55 of 98 neo-epitopes get a negative ExpectedTpm (a pan-cancer median of 0.382 is written as -0.618;
    UBR3's 14.56 as 13.56). Private pair: 18 of 44.
  - Peptides (HG008): EffectiveTpm changes on every row, ExpLikelihood on 624 rows (most are below the printed
    precision), ExpLikelihoodRank on 26,599.
  - Corrected: a sum starts at 0 when its first value is added; without a value it stays -1.
- **NE7 (scorer, no effect).** TpmCalculator tests `x.NeIds.contains(neoData1)`: a Set<Integer> against a
  NeoEpitopeData, never true. So a neo-epitope already in a shared upstream group starts a second, smaller group.
  - Because matching transcripts is an equivalence (the lists hold distinct names), the later group is a subset of
    the earlier one, and the lookup always finds the earlier.
  - Corrected (the id is tested). The unit test `ne7_shared_upstream_groups` gives equal TPMs in both modes.
- **NE8 (scorer).** ScoreConfig.getScoringFilename prefixes `-score_file_dir` to a score-file option's value, though
  ConfigBuilder has already checked that value as an existing path as given.
  - So a full path, or any path not also inside the directory, is not found. `sc_dir_and_file`: the recognition
    file is not found and Java exits 1.
  - Corrected: the value is used as given.
- **NE9 (scorer).** `-peptide_lengths` ("Peptide length min-max, separated by '-', eg 8-12") is registered but never
  read: the peptides are always 8 to 12 long.
  - Corrected: the range is applied. An invalid value is an error; lengths outside 8-12 get no binding likelihood,
    as in Java.
  - The pipeline does not pass the option.
- **NE10 (scorer).** With a sample id file and `-threads` > 1, Java splits the samples between tasks that write to the
  shared files under a lock, so the rows of the samples interleave in thread order.
  - Shown by `sc_cohort_threads2`: Java's files differ from those at 1 thread (the rows are the same).
  - The port writes the samples in order at any thread count, which is Java's 1-thread output in both modes (the
    defined reference; no switch).
- **NE11 (scorer).** `RnaSomaticVcf = getValue(RNA_SOMATIC_VCF, sampleDataDir)`: with `-sample_data_dir` and no
  `-rna_somatic_vcf`, the RNA VCF is the sample directory.
  - Java tries to read the directory as a VCF and gets a null variant list. setMutationRnaSupport then throws on the
    first point mutation (exit 1, `sc_sample_dir_only`).
  - So every DNA-only run that uses `-sample_data_dir` fails.
  - Corrected: the RNA VCF is read only when given. HasRna already requires the option.
- **NE12 (scorer, kept: intent unclear).** Samples from a sample id file get a null RNA sample id.
  - With HasRna true, the ISOFOX files are looked up as `null.isf.transcript_data.tsv` and `null.isf.neoepitope.tsv`
    (`sc_rna_sample_file` provides such files), and no RNA genotype is found in the VCF.
  - The unused constant RNA_SAMPLE_ID_SUFFIX ("_RNA") suggests the intended id, but that is not certain, so Java's
    behaviour is kept in both modes.
- **NE13 (scorer).** BindScorer.calcScoreData passes the score rank to getBindingLikelihood even when it is
  INVALID_SCORE (-1, no score distribution for the allele and length). -1 is below the first rank bucket, so Java
  returns the highest likelihood of that allele and length.
  - Shown by `sc_no_rank_dist`, the score files without A0201's distribution: Java gives all 10,724 A0201 rows rank
    -1 and the top likelihood, and 8,409 of them a likelihood rank below 0.02 (171 with the distribution). They look
    like strong binders.
  - Corrected: an invalid rank gives an invalid likelihood, as an unknown allele or length does. The rows then get
    likelihood -1 and likelihood rank 1, Java's own handling of an invalid likelihood.
  - The bundled files have a distribution for every allele with a matrix, so this is latent in the pipeline.
- **NE14 (finder).** NeoCommon.IMMUNE_TRANSCRIPT_PREFIXES is {TR_J, TR_V, TR_D, IG_J, IG_D, IG_J}: IG_J twice and no
  IG_V. NeoSampleTask's "ignore IG and TR transcripts" therefore keeps the IG V-segment transcripts.
  - The Ensembl cache has 143 coding IG_V genes.
  - Shown by `fin_ig_v`: Java writes a neo-epitope for a synthetic missense in IGKV1D-17.
  - Corrected: IG_V in place of the second IG_J. No comparison sample has such a mutation.

### Kept as Java (oddities with no effect worth a finding, or choices)

- BindData.setTPM tests the current TPM (0 at first) rather than the new one, so negative TPMs are stored. They give
  the lowest expression bucket, as 0 would. `hasTPM()` is always true.
- AlleleDisrupted and AlleleCN come from the first LILAC row of an allele. For a homozygous allele, the second copy's
  copy number and variants are not seen.
- A peptide shared by neo-epitopes of a variant is written once, with the first such neo-epitope's type, info and
  gene, and every id in NeIds. A peptide found in one of them's upstream wild-type protein is dropped from that one
  only.
- With RNA, ExpectedTpm uses the local copy number; without, the sample ploidy (setCalculatedTpmValues). This is a
  modelling choice.
- `-write_types PEPTIDE` and `SAMPLE_SUMMARY` are accepted and write nothing. `-write_peptide_scores`,
  `-write_summary_data` and `-validation_data_file` belong to BindScorer's own main class.
- Without `-ensembl_data_dir`, or without the amino acid file, Java ignores the failed load and drops no wild-type
  peptides. The port does the same.
- `-rna_somatic_vcf` is described as a "Directory for Purple somatic variant RNA-appended files", but it is a file
  path (with `*` for the sample id).
- No sample (`-sample` and `-sample_id_file` absent): Java creates `neo_cohort.*.tsv` in the output directory and
  exits 0 leaving them empty. The port does the same.
- An output directory that does not exist: Java logs an error and exits 0 without outputs. The port does the same.
- The port reads only the TPM-medians columns the samples need ("All" and their cancer types), where Java parses
  every column. They differ only for a malformed value in another column.
- EpitopeUtils.generatePeptides' `(int)round(fragmentCount)` of an integer, and the Poisson range cache from the
  bundled `expected_fragment_ranges.csv` (0 to 9 fragments), are as Java.

## Tests

- **neo-rs:** `CARGO_BUILD_JOBS=2 cargo test --release`, 32 tests, all passing.
  - NEO's Java tests of the ported classes, translated by a script and finished by
    hand; every assertion kept except AminoAcidRna's DNA/RNA swaps (not ported):
    - PmNeoEpitopeTest: SNV, MNV, deletion, insertion and stop-lost;
    - SvNeoEpitopeTest: fusions and the non-disruptive SV test;
    - NeoEpitopeUtilsTest: DNA routines, upstream and downstream coding bases, skipped acceptors and donors, peptide
      combinations;
    - TpmCalcsTest: effective TPMs, complex shared peptides, Poisson ranges.
  - The test helpers follow hmf-common's GeneTestUtils.createTransExons and MockRefGenome.
  - Tests of my own:
    - NE6, NE7, NE8, NE9, NE11, NE13 and NE14 in both modes;
    - the score config, file names and sample id files;
    - score ranks in ascending and descending distributions;
    - binding likelihoods with the two-digit and gene fall-backs;
    - expression likelihoods; position and flank weights;
    - recognition similarity, and that the per-peptide similarities equal Java's per-allele calls;
    - a peptide's neo-epitope id order (Guava's two-bucket HashSet);
    - the readers (neo-epitopes, LILAC, ISOFOX, the RNA VCF);
    - TPM normalisation, effective TPM and the copy number factor.
- **Not ported:** BindModelTest and BindUtilsTest, which test the training classes (BindTrainer, PosWeightModel,
  NoiseModel, RandomDistributionTask), not used by NeoScorer; MissensePeptideTest (the missense main class).
- **hmf-common-rs:** `cargo test --release`, 112 tests, all passing.
  - Session 1's `neo.rs`, `hla.rs` and `rna.rs` tests are run here for the first time.
  - Added: `ensembl::tests::transcript_amino_acids` and `poisson::tests::poisson_noise_values`.

## Speed and memory

At one thread the port is 2.4 to 2.8 times faster than Java for the finder and about 10 times for the scorer, in an
eighth to a tenth of the memory (table above).

- **The finder** (HG008, 0.72 s): mostly the Ensembl cache (session 1).
- **The scorer** (HG008, 0.80 s at 1 thread, 0.64 s at 2; `-log_debug` gives the phase times):
  - reference data 0.20 s: the TPM medians, 215k transcripts from gzip, 0.11 s; the Ensembl transcript amino acids
    0.08 s;
  - the binding reference files 0.06 s;
  - TPM calculator, 10,789 peptides: 0.004 s;
  - scoring 128,688 allele-peptides: 0.38 s at 1 thread, 0.22 s at 2. Most of it is the recognition similarity
    (BLOSUM distances to every allele's immunogenic peptides of the same length), found once per peptide rather than
    once per allele as in Java;
  - writing 64k rows: 0.15 s.
- **Memory** (175 MB) is mostly the medians and the transcript amino acids. Java: 1.4 GB.
- The pipeline gives the scorer no `-threads`, so it runs at 1 thread there; the finder uses threads only for several
  samples.

## Shared code: additions to hmf-common-rs

Everything is an addition; no existing function changed behaviour.

- **Session 1:**
  - `neo.rs` (new): NeoEpitopeType; NeoEpitopeFile (the neo_data.tsv header, line and reader); NeoEpitopeFusion
    (LINX's neoepitope reader); RnaNeoEpitope (ISOFOX's neoepitope reader); pointMutationInfo, fusionInfo,
    extractTranscriptNames;
  - `hla.rs` (new): the LilacAllele reader (`<sample>.lilac.tsv`, with the .csv fall-back);
  - `rna.rs` (new): the ISOFOX TranscriptExpressionFile reader and RnaExpressionMatrix (the cohort matrix);
  - `ensembl.rs` (added): TranscriptAminoAcids, `ensembl_amino_acid_data_file`, `load_transcript_amino_acid_data`;
  - `poisson.rs` (added): PoissonCalcs (`calc_poisson_noise_value`, `calc_poisson_noise_value_with`,
    `find_upper_start_value`);
  - `lib.rs`: the three new modules.
- **Session 2:** tests only (`ensembl.rs`, `poisson.rs`, above).

## Open points

- **NE12:** the RNA sample id of sample-id-file samples. If Hartwig confirms the intent (`<id>_RNA`, or the tumor
  id), it is a one-line correction.
- **NE4, NE13, NE14:** reached only by synthetic inputs or unit tests. None of the samples has a non-coding upstream
  fusion transcript, an allele without a score distribution, or an IG V-segment mutation.
- **In the pipeline:** the port has not yet run as NEO_FINDER / NEO_SCORER in an oncoanalyser run (the
  `--rust_ports_bin` profile needs `neo-rs` added; the wrapper's main-class argument is supported).
- **Not ported:** NEO's other main classes (BindTrainer, BindScorer, MissensePeptideScorer, the `utils` programs).
