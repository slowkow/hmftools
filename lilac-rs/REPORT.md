# LILAC in Rust: record of the port

A record of how `lilac-rs` was checked against Java LILAC 2.1 (this repository's `lilac/target/lilac-2.1-jar-with-dependencies.jar`),
the faults found in the Java code, and timings. Options and outputs are in [README.md](README.md).

## 1. How it was checked

Java LILAC 2.1 and the port were run on the same inputs with the same options, and every output file was
compared byte for byte. Java ran on the local JDK 23 (aarch64) with `-Xmx8g -XX:+UnlockExperimentalVMOptions
-XX:hashCode=2`: with constant identity hashes, the HashSets of enum constants that set the order of the gene
groups (LI1) and of the QC statuses (LI8) iterate in insertion order, which is the defined reference the port
reproduces in both modes. The port ran in Java mode (`LILAC_JAVA_COMPAT=1`).

Inputs: a private whole-genome tumor/normal pair (ids TUMOR / NORMAL; LILAC reads only the HLA genes'
coding regions), renamed copies of the PURPLE outputs of an earlier pipeline run on a second private pair, another
person (gene copy numbers, somatic VCF; they exercise the code paths), the public GIAB HG008 HLA-region slices (chr6:29-34 Mb, tumor cell line 71x, normal 30x) with that
run's PURPLE outputs, and a synthetic somatic VCF for HG008 (non-reference SNVs of the HG008 tumor in
HLA-A exons taken from a pileup at run time, one marked SYNONYMOUS, one filtered, plus a frameshift deletion:
made by a comparison script kept outside this repository), since neither real VCF has a variant in the HLA coding regions. The
pipeline's options: `-genes ALL -freq_score_penalty 0.0009 -sequencing_type ILLUMINA -threads 2`, or for
panels `-genes MHC_CLASS_1 -freq_score_penalty 0.0018 -targeted_panel`.

| case | options besides the common ones | Java time, peak RSS | Rust time, peak RSS | result |
|---|---|---|---|---|
| wgs_tn | WGS T/N, -purple_dir | 29.0 s, 5.25 GB | 0.72 s, 0.40 GB | identical |
| wgs_tn_t1 | the same, -threads 1 | 29.6 s, 5.25 GB | 1.35 s, 0.32 GB | identical (Java 1 vs 2 threads: identical) |
| wgs_tn_nopurple | WGS T/N, no PURPLE | 28.6 s, 5.32 GB | 0.66 s, 0.38 GB | identical |
| wgs_go | germline only (normal) | 26.9 s, 5.36 GB | 0.56 s, 0.34 GB | identical |
| wgs_to | tumor only, -purple_dir | 27.2 s, 5.14 GB | 0.55 s, 0.38 GB | identical |
| wgs_c1 | -genes MHC_CLASS_1 | 12.2 s, 4.68 GB | 0.69 s, 0.28 GB | identical |
| wgs_panel | panel settings (class I, -targeted_panel, 0.0018) | 12.6 s, 4.63 GB | 0.68 s, 0.28 GB | identical |
| wgs_all | -write_types ALL -log_debug (all detail files) | 30.4 s, 5.24 GB | 0.71 s, 0.40 GB | identical (20 files) |
| hg008_tn | HG008 T/N, -purple_dir | 32.0 s, 5.30 GB | 0.96 s, 0.29 GB | identical |
| hg008_tn_t1 | the same, -threads 1 | 38.0 s, 5.31 GB | 1.11 s, 0.23 GB | identical (Java 1 vs 2 threads: identical) |
| hg008_go | germline only | 31.0 s, 5.38 GB | 0.57 s, 0.26 GB | identical |
| hg008_to | tumor only | 32.6 s, 5.17 GB | 0.86 s, 0.30 GB | identical |
| hg008_c1 | -genes MHC_CLASS_1 | 15.9 s, 4.92 GB | 0.64 s, 0.23 GB | identical |
| hg008_panel | panel settings | 16.4 s, 5.05 GB | 0.65 s, 0.23 GB | identical |
| hg008_rna | -rna_bam (the tumor slice as a stand-in) | 33.1 s, 5.19 GB | 0.71 s, 0.29 GB | identical |
| hg008_all | -write_types ALL | 33.0 s, 5.40 GB | 0.67 s, 0.28 GB | identical (21 files) |
| hg008_somatic | the synthetic somatic VCF + gene copy numbers | 32.3 s, 5.38 GB | 0.90 s, 0.29 GB | identical (3 variants assigned, 1 not: WARN_UNMATCHED_SOMATIC_VARIANT) |
| hg008_sbx | class I, -sequencing_type SBX (1,652,400 complexes) | 20.6 s, 5.27 GB | 0.86 s, 0.25 GB | identical |
| hg008_ultima | class I, -sequencing_type ULTIMA (no UQ tags: no low-qual bases; 5,103,000 complexes) | 51.1 s, 6.12 GB | 3.51 s, 0.36 GB | identical |
| hg008_nocap | class I, -max_ref_fragments 0 (no down-sampling) | 15.3 s, 5.04 GB | 0.62 s, 0.23 GB | identical |
| hg008_cap300 | class I, -max_ref_fragments 300 (every n-th fragment) | 14.7 s, 4.95 GB | 1.40 s, 0.24 GB | identical |
| hg008_drb | -genes HLA_DRB1,HLA_DRB | 9.4 s, exit 1 | 0.14 s, exit 1 | same failure (LI23), empty files |
| hg008_c1_t0 | class I, -threads 0 | 10.7 s, 4.73 GB, exit 1 | 0.57 s, exit 1 | same failure (LI7), empty files |

(Java peaks reflect -Xmx8g. The machine was shared with other jobs (load 6-11); Rust times vary by up to
1.5x with load: by instructions retired, -threads 1 and 2 cost the same, 16.9 G for wgs_tn.)

In the default (corrected) mode the typing outputs `lilac.tsv` and `lilac.qc.tsv` equal Java's in every case;
`candidates.coverage.tsv` differs only in its header line (LI2), and the REF_COUNTS files of HLA-B and HLA-C
in the raw depths at exon-boundary loci (LI3).

Hash and thread dependence of Java: with the JVM's default identity hashes, the order of the gene group rows
changes from run to run (LI1; two runs here: DPB1, DQA1, CLASS_1, DQB1, DPA1; the pipeline's run: DPB1, DPA1,
CLASS_1, DQB1, DQA1; -XX:hashCode=2: CLASS_1, DQB1, DPA1, DPB1, DQA1), the rows themselves are identical. With
10,000 complexes or more the result can depend on the thread count (LI6): on HG008 (982,800 complexes)
Java's two 2-thread tasks cull 480,492 and 481,507 complexes, one thread culls 962,002, so the ranked sets
differ, but the kept complexes and the outputs are identical. Everything else is deterministic: string
hashes (allele names, read names) fix the HashMap orders, and the gene tasks' results are joined in order.

## 2. Findings

Each correction is behind `compat::java("LIn")`; `-java_compat` restores Java's behaviour. "Reached" says
whether the case set above reaches the code (a comparison script counts the `LIn:` debug traces).

| id | what Java does | why it is wrong | evidence, effect | port |
|---|---|---|---|---|
| LI1 | solves and writes the gene groups in the order of a `HashSet<GeneSelector>` (identity hashes) | the row order of all three output files changes between runs | two default-JVM runs and the pipeline's run differ in order, rows equal | always the -XX:hashCode=2 order (as given, ALL in enum order) |
| LI2 | opens `candidates.coverage.tsv` with the header of the fragment file (`Complex FragmentId FragmentCoords Type FullAlleles WildAlleles`); `HlaComplexFile.header()`, which matches the rows, is never used | the header names 6 columns for 14-column rows | every run | corrected: the matching header |
| LI3 | the phasing fragments of each gene are enriched (splice bases added) in place on the "original" copies, which the comment says stay unchanged; later genes' phasing and raw counts then see them | gene B's and C's phasing depends on gene A's | WGS: the raw depth of an HLA-B exon-boundary locus 75 instead of 18 (REF_COUNTS files); typing unchanged here; can change the low-depth filter (10) in thin data | corrected: the enrichment works on copies |
| LI4 | `getReferenceNucleotides` lists the nucleotides of the loci with counts by position and is read by locus | a locus without counts shifts every later locus by one | not reached (no gaps) | corrected: indexed by locus |
| LI5 | `Candidates.filterCandidates` deletes low-depth loci from the phased evidence's own loci list, but the evidence sequences keep them | the haplotype QC then compares misaligned sequences: every sequence looks unmatched and haplotypes are built from the wrong amino acids; the amino acid QC skips the loci of these haplotypes | reached in every run: HG008 5 unmatched haplotypes instead of 0, WGS 9 instead of 2; none had 3 supporting fragments, so the QC columns agree here | corrected: the filter works on copies |
| LI6 | with -threads > 1 and 10,000 complexes or more, deals the complexes out to the threads in turn; each thread culls by its own running maximum, and the shares are joined thread by thread | the complexes kept, and their order into a non-transitive sort (LI19), depend on the thread count | HG008: 962,002 culled at 1 thread, 961,999 at 2 (above); same result here | corrected: single-thread order and culling at any thread count (the totals are still computed in parallel) |
| LI7 | -threads 0 (described as "0 or 1 not multi-threaded") builds `Executors.newFixedThreadPool(0)` once a group has 10,000 complexes: IllegalArgumentException, exit 1, empty files | crash | hg008_c1 -threads 0 | corrected: 0 taken as 1 |
| LI8 | the QC `Status` column joins a `HashSet<LilacQCStatus>` | the order of several statuses changes between runs | (only one status in these cases) | always insertion order |
| LI9 | sorts the somatic coding counts by `HlaAllele.compareTo` (protein number as a number) and pairs the i-th with the i-th allele of the coverage, sorted by name | e.g. A*02:11 and A*02:101 sort differently, so the counts go to the other allele | not reached | corrected: each allele's own counts |
| LI10 | `writeVariant` returns when the variant context is not null, which it never is | the `<sample>.lilac.somatic.vcf.gz` the code builds is never written | Java never writes it | not corrected: no consumer in the pipeline (open point) |
| LI11 | `PhasedEvidence.equals` returns false when any locus matches and true when all differ (with equal hash: List.hashCode collisions such as [10,50] and [11,19]) | two different evidences can be taken as one in the HashSet of the merge step | not reached | corrected: identity, as the code uses it |
| LI12 | `expandToSixAlleles` gives a homozygous allele's two copies half the unique and shared coverage each, but the wildcard coverage only to the first | the second copy's wildcard (and total) coverage is lost | not reached (no wildcard coverage) | corrected: the remainder |
| LI13 | the debug file `HLA_Y.fragments.tsv` is never closed | only whole 8192-character buffers reach the disk | debug mode with HLA-Y above threshold only | corrected: written whole |
| LI14 | `HlaAllele.compareTo` returns 0 for alleles whose names have the same String hash | four pairs in the resources (A*32:11Q/A*32:122, C*03:22Q/C*03:232, C*16:16Q/C*16:172, A*02:01:14Q/A*02:01:152): a pair can be ordered both ways (a duplicate complex) and ties break by input order | not reached | corrected: compared by fields |
| LI15 | `ExpectedAlleles`: index i holds the count for locus i+1 but is read for locus i | an apparent off-by-one at the boundaries where genes become distinguishable | intent unclear | kept |
| LI16 | `HLA_Y.coverage.tsv` lists the sources in `HashMap<FragmentSource>` order | identity-hash order (REF_COUNTS, HLA-Y above threshold) | | insertion order |
| LI17 | REF_COUNTS: `candidates.aminoacids.txt` and `.nucleotides.txt` have no gene group in their names | each group overwrites the previous: only the last group's counts remain | wgs_all, hg008_all | kept |
| LI18 | `createAlignmentFragments` builds one fragment per alignment block (block end one past) but each from the whole read, then merges them | the blocks only decide whether there is a fragment | no effect | kept |
| LI19 | `ComplexCoverageSorter` treats scores within 0.0001 as equal | not transitive: the order depends on the input order, and TimSort may in principle throw | | TimSort reproduced |
| LI20 | the result depends on the JVM: `HashSet(Collection)` is sized differently since JDK 19, and `Math.log10` is an intrinsic on x86-64 | complex order and score ties may differ between JDKs and platforms | | JDK 19+ sizes, fdlibm log10 (= the aarch64 JVM: 2,290 of 2,290 frequency logs equal; Rust's own log10 differs on 26) |
| LI21 | option texts: `-actual_alleles` / `-restricted_alleles` "comma separated" but split at `;`; `-write_types` offers READS, which is rejected; `MedianBaseQuality` is the 25th percentile | misleading | | kept |
| LI22 | `AlleleCoverage.increment` finds the unique-count key by identity | another object of the same allele would replace the count | no effect (cached alleles) | kept |
| LI23 | `-genes HLA_DRB1` (not in ALL) on HG008: all candidates are unsupported wildcard alleles, 0 complexes, "failed to calculate complex coverage", exit 1 | the DRB groups are not usable on this sample | | the same failure |
| LI24 | the trimmed HLA-DRB1 sequences keep the wildcard flag of the untrimmed sequence | an allele can count as having wildcards it no longer has | DRB1 only | kept |
| LI25 | with down-sampled reference fragments, the unsupported wildcard alleles are removed only from the sampled fragments; the others keep them, and fragments left with only those count as CANDIDATE not WILD_ONLY | inconsistent, both count as unmatched | no QC effect | kept |

## 3. Java 2.1 and the pipeline's Java 2.0

HG008: the pipeline's LILAC 2.0 output (`data/public/hg008/pipeline/lilac/`) and Java 2.1 on the HLA slices
have the same rows in all three files (A*02:01, A*03:01, B*07:02, B*15:01, C*03:04, C*07:02; DQB1*03:02,
DQB1*06:02; DPA1*01:03 x2; DPB1*04:01 x2; DQA1*01:02, DQA1*03:01), in another gene group order (LI1). The
tumor's copy numbers of 0.98 and 0.03 per gene are the cell line's loss of one chr6p copy with purity 1.0.

Private pairs: the earlier pipeline run behind the PURPLE files is of another person. Its Java 2.0
types share exactly one allele per gene with those of the BAMs used here, whose tumor and normal type
identically to each other. That is how the two pairs were found to be different people.

## 4. Speed

Java spends most of its time reading the two sequence files again for every gene group (about 4.5 s per
group); the port reads them once (0.08 s) and parses each group's lines. The rest, for the whole-genome pair
at 2 threads: BAM reading and fragment building 0.07 s per BAM, phasing 0.07 s, the fragment-allele mapping
0.2 s (parallel), the complexes 0.08 s (HG008 with 982,800 complexes: 0.15 s; Java 5.6 s), the tumor side
0.1 s; gene groups run in parallel. With -threads 1 the whole-genome pair takes 1.3 s, with 2 threads 0.67 s.
