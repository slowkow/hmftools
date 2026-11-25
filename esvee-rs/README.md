# esvee-rs

A Rust port of [ESVEE](../esvee/README.md) 2.1, the structural variant caller. One binary holds the four
steps of the Java tool (prep, assembly, depth annotation, calling); each step takes the command line of
the Java application of that name and writes the same files.

By default the port corrects 63 faults found in the Java code while porting, so its output differs
from Java's where those faults matter. With `-java_compat` it behaves as Java does and its output is
identical to Java's; that is how the port was verified. See "Corrections and `-java_compat`" below.

Status and measurements are in [PLAN.md](PLAN.md). What has been compared with Java, and what Java does
that is worth knowing, is in [REPORT.md](REPORT.md).

## Build

```
cargo build --release
```

The binary is `target/release/esvee-rs`.

- A C compiler is needed. `build.rs` compiles the BWA-MEM sources in `vendor/bwa` (through the `cc`
  crate) and links them into the binary; the assembly step aligns with them. zlib must be available to
  the linker.
- htslib is linked statically (through `rust-htslib`), so samtools is not needed at run time.
- The default build uses the mimalloc allocator; `cargo build --release --no-default-features` uses the
  system one.

All runs recorded in REPORT.md were made on macOS on Apple silicon. On ARM the SSE2 code of BWA is
built through `sse2neon.h`. No other platform has been tried yet.

## Run

```
esvee-rs <prep|assembly|depth|caller> <the Java step's options>
```

Options have the names and defaults of the Java step (`PrepApplication`, `AssemblyApplication`,
`DepthAnnotator`, `CallerApplication`). The four steps in order, for a tumor/normal pair:

```
esvee-rs prep \
    -sample TUMOR,NORMAL -bam_file tumor.redux.bam,normal.redux.bam \
    -ref_genome ref.fasta -ref_genome_version 38 \
    -known_hotspot_file known_fusions.38.bedpe \
    -output_dir out/ -threads 8

esvee-rs assembly \
    -tumor TUMOR -reference NORMAL \
    -tumor_bam out/TUMOR.esvee.prep.bam -reference_bam out/NORMAL.esvee.prep.bam \
    -junction_file out/TUMOR.esvee.prep.junction.tsv \
    -ref_genome ref.fasta -ref_genome_version 38 \
    -output_dir out/ -threads 8

esvee-rs depth \
    -sample TUMOR,NORMAL -bam_file tumor.redux.bam,normal.redux.bam \
    -input_vcf out/TUMOR.esvee.raw.vcf.gz \
    -ref_genome ref.fasta -ref_genome_version 38 \
    -unmap_regions unmap_regions.38.tsv \
    -output_dir out/ -threads 8

esvee-rs caller \
    -tumor TUMOR -reference NORMAL \
    -input_vcf out/TUMOR.esvee.ref_depth.vcf.gz \
    -ref_genome_version 38 \
    -pon_sgl_file sgl_pon.38.bed.gz -pon_sv_file sv_pon.38.bedpe.gz \
    -known_hotspot_file known_fusions.38.bedpe \
    -repeat_mask_file repeat_mask_data.38.fa.gz \
    -output_dir out/ -esvee_prep_dir out/
```

The assembly step needs the BWA index image of the reference genome beside the FASTA (`<fasta>.img`),
as the Java step does. It is the same file; nothing has to be rebuilt for the port.

Outputs, with the Java names: `<sample>.esvee.prep.bam` (+ `.bai`) per sample and
`<tumor>.esvee.prep.junction.tsv`, `.prep.fragment_length.tsv`, `.prep.disc_stats.tsv`;
`<tumor>.esvee.assembly.tsv`, `.breakend.tsv`, `.raw.vcf.gz` and the optional assembly TSVs of
`-write_types`; `<tumor>.esvee.ref_depth.vcf.gz`; `<tumor>.esvee.unfiltered.vcf.gz`, `.somatic.vcf.gz`
and `.germline.vcf.gz`, each VCF with a tabix index.

### All four steps in one call

Without a step name, esvee-rs does what the Java `EsveeApplication` does (the class the bioconda
`esvee` wrapper runs, and so the command nf-core/oncoanalyser runs): it takes the options of all four
steps on one command line and runs prep, assembly, depth and caller in turn.

```
esvee-rs \
    -tumor TUMOR -tumor_bam tumor.redux.bam -reference NORMAL -reference_bam normal.redux.bam \
    -esvee_prep_dir out/ -ref_genome ref.fasta -ref_genome_version 38 \
    -known_hotspot_file known_fusions.38.bedpe -pon_sgl_file sgl_pon.38.bed.gz -pon_sv_file sv_pon.38.bedpe.gz \
    -repeat_mask_file repeat_mask_data.38.fa.gz -unmap_regions unmap_regions.38.tsv \
    -sequencing_type ILLUMINA -bamtool $(which sambamba) \
    -write_types 'PREP_JUNCTION;PREP_BAM;FRAGMENT_LENGTH_DIST;JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF' \
    -threads 8 -output_dir out/
```

Leave out `-reference` and `-reference_bam` for a tumor-only run. `-decoy_genome`, `-saga_fasta` and
`-target_regions_bed` are taken as by the Java steps.

- **Options.** The set is the union of the four steps' options, with Java's names, types and defaults
  (`esvee-rs -help` lists them, in the format of Java's ConfigBuilder). An option a step registers goes
  to that step; an unknown option is an error, as in Java. Path options are checked when the run
  starts, before prep has written anything, so `-esvee_prep_dir` must exist beforehand (oncoanalyser
  runs `mkdir -p esvee/`). Java also needs `-output_dir` to exist; the port creates it.
- **How each step finds its inputs** is that of the Java configs when they share one command line:
  prep and depth read `-sample` and `-bam_file`, or else the tumor and then the reference ids and BAMs;
  assembly ignores `-tumor_bam` and `-reference_bam` and reads the prep BAMs and junction file from
  `-esvee_prep_dir` (or `-output_dir`); depth reads `<output_dir><tumor>.esvee.raw.vcf.gz`; the caller
  reads `<output_dir><tumor>.esvee.ref_depth.vcf.gz` and the prep files from `-esvee_prep_dir`. An
  `-input_vcf` given on the command line is read by both depth and the caller, as in Java.
- **`-write_types`** is read by prep for the prep types and by assembly for the assembly types. Prep
  writes `disc_stats.tsv` in any case.
- **`-stages PREP;ASSEMBLY;REF_DEPTH;CALLING`** runs only the stages listed (in that order, whatever the
  order given), the later ones reading the files of earlier runs.
- **`-version`** prints `Esvee version 2.1` and exits. (Java does this only when `-version` is the only
  argument.)
- `-java_compat` and `-java_compat_except` apply to all four steps.
- **`-bamtool`**: the prep BAMs are sorted and indexed in the process with either tool named (see
  "Differences from the Java tool"). Java's `sambamba sort` gives the same records in the same order as
  `samtools sort`; only the `@PG` lines of the header differ.

The output equals that of the four steps run one at a time with the matching options, and with
`-java_compat` that of Java's `EsveeApplication` (REPORT.md, section 2.5). Run Java's reference for such a
comparison as for the steps (one thread, `-XX:hashCode=2`); in one JVM the option also applies to prep,
which then takes about ten times as long (its identity-keyed sets degrade to lists) with the same output;
`-stages PREP` without the option and then `-stages 'ASSEMBLY;REF_DEPTH;CALLING'` with it avoids that.

## What is ported

| Step | Ported | Not ported |
|---|---|---|
| prep | everything: tumor/normal, tumor only, normal only, BAM and CRAM input, Illumina, SBX and Ultima, the cache BAM and `-no_cache_bam` modes, `-write_types`, SAGA location matching (`-saga_fasta`) | nothing known |
| assembly | everything: the TSVs, the raw VCF and the visualiser files; tumor/normal, tumor only, germline only, Illumina, SBX and Ultima, SAGA, `-decoy_genome`, all `-write_types`, `-perf_debug` and `-perf_log_time` log lines | nothing known; the visualiser (`-run_assembly_vis` with `-specific_junctions`) writes the same `vis/*.html` as Java on the machine it was checked on (text positions come from that machine's font metrics) |
| depth | everything, BAM and CRAM (3.0 and 3.1) | htsjdk's record validation behind `-bam_validation STRICT` and `LENIENT`: the option is accepted, records are not validated |
| caller | everything: tumor/normal, tumor only, germline only, PON and artefact PON, hotspots, repeat mask, target regions, `-write_breakend_tsv`, `-manual_ref_depth` | nothing known |

What has and has not been tested for each is in REPORT.md, sections 2 and 6. In short: with `-java_compat` all four steps
give Java's output on a small Illumina tumor/normal input and across the option cases, and on an
Illumina whole-genome tumor/normal pair every output from assembly on is identical to Java's (prep:
the same records, with one pair of tied records per BAM in the other order, as Java itself does from
run to run); reference genome 37 has not been run on real data in any step.

## Corrections and `-java_compat`

The faults found in Java ESVEE are numbered in REPORT.md, section 4. For 63 of them the port takes a
corrected path by default; section 4.5 says what each correction does and what it changes.

- `-java_compat` on any step (or `ESVEE_JAVA_COMPAT=1` in the environment) gives Java's behaviour for
  every finding. Use it to compare with Java, or to continue an analysis made with Java ESVEE.
- `-java_compat_except A12,P10` gives Java's behaviour for every finding but those listed.
- The corrections that change the most on real data are P10 (short supplementary reads in prep) and
  A12 (adapter trimming in assembly). For both the intent of the Java code is not certain; the readings
  taken are described in REPORT.md and are open to correction by the ESVEE authors.
- The corrected output cannot be checked against Java. It is checked by unit tests of the corrected
  paths, by being the same at any thread count, and by measuring what each correction changes.

## Base-quality bands (options of the port, not of Java)

Java ESVEE fixes two quality bands: a base at Q26 or below is low quality (prep drops a read with more
than half its bases low, and marks a soft clip with 5 or more low bases; assembly lets low bases
mismatch), and an Illumina base at Q30 or above is high quality. On data whose base qualities were
rewritten by recalibration and stop below Q30, most of a tumor's bases can be low, and its junctions
then lose their support. `-low_base_qual_threshold N` and `-high_base_qual_threshold N` move the bands
for prep and assembly (defaults 26 and 30; the low band must end below the high one); in the one-call
mode they are given to both steps, and `-help` lists them after Java's options. Not given, the output is
unchanged, with or without `-java_compat`; given, the log says so. Hartwig tuned the later filters on
data with the usual qualities, so calls made with moved bands need checking.

## What "the same" means

All of this section is about the port run with `-java_compat`.


- **BAM:** the same records in the same order with the same fields and attributes, and the same header
  apart from `@PG` lines, compared after decompression.
- **TSV:** byte-identical.
- **VCF:** byte-identical after decompression, header included.
- **Compressed bytes and index files** (`.bai`, `.tbi`) differ from Java's, since the compressor and
  the index writer differ. The indexes were checked by querying through them.
- **Thread count:** the port writes the same output at any thread count, and from run to run. (Until
  2026-11-24 this failed on reads from sequencing runs with read-name prefixes of different lengths, as
  when tumor and normal were sequenced separately: assembly cut read names as Java's `ReadIdTrimmer`
  does, finding A5, which is now corrected by default. With `-java_compat` that fault is kept, and on
  such data assembly output then depends on thread timing, as Java's does.) Java does not:
  - Java's prep output depends on the thread count (the partitions share a cache of spanning reads).
    The port gives Java's one-thread result.
  - Java's assembly output depends on the thread count and also changes from run to run at one thread,
    because sets of objects are iterated in identity-hash order. The port gives the result of Java at
    one thread with every identity hash made equal (`-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`).
  - Java's depth and caller VCFs do not depend on the thread count, and the port matches them.

  To compare the port with Java, run the port with `-java_compat` and Java with `-threads 1` and, for
  the assembly step, with that JVM option. REPORT.md, section 3, has the detail.

## Differences from the Java tool

- **Prep sorts in the process.** Java writes an unsorted BAM and calls `samtools sort` and
  `samtools index` through `-bamtool`. The port sorts by the same key (ties keep the order in which
  the records were written) and builds the index itself. `-bamtool` only switches sorting on, as in
  Java: without it the unsorted BAM is left. The sorted BAM's header has `SO:coordinate` and no `@PG`
  line for samtools. Sort runs of 768 MB spill to temporary BAMs; the environment variable
  `ESVEE_PREP_SORT_RUN_MB` changes the size.
- **`final_breakend.tsv` (caller, `-write_breakend_tsv`):** Java writes the filters of a row in an
  order that changes from run to run; the port writes them in the order of the `FilterType` enum.
- **Decoy rows (assembly, `-decoy_genome`):** always in junction group order, which is Java's order at
  one thread; Java at more threads writes them as tasks finish.
- **CRAM 3.1** is read (depth step tested); the htsjdk in the Java tool reads 3.0 only.
- **Failures.** Where Java stops with an uncaught exception the port logs one error line and exits
  with status 1. The wording differs, and in some cases so does what is left on disk (for example a
  complete header-only VCF where Java leaves a truncated file).
- **Faults in Java are corrected by default** and reproduced with `-java_compat` (above). Each place
  is marked in the code with the finding's number (`compat::java("A12")`).
- **A missing SAGA index image** stops the assembly step at start-up; Java only fails when the first
  matcher is created.
- **Log lines** have the Java text. Timings in them differ, and with `-perf_log_time` fewer items
  exceed a given threshold. The caller prints the lines about its reference data together once the
  data is loaded.
- **Memory** is lower in every step on the small input (REPORT.md, section 5). On a whole genome the
  assembly step peaks at about 7 GB, plus the mapped BWA index image of 5-6 GB.
- **`ESVEE_MEM_STATS=1`** makes the assembly step log its memory use and live object counts at each
  stage.

## Tests

`cargo test --release` runs the unit tests: the Java unit tests of the four packages ported under their
names, tests of the emulated Java collections against values taken from a JVM, and tests of the port's
own code. The caller's tests include end-to-end runs on synthetic inputs in `src/caller/testdata`
against the files Java wrote for them.

The comparisons with Java on real data are driven by scripts kept outside the repository, because the
outputs they store are derived from patient data.

## Licence of the vendored BWA

`vendor/bwa` holds the alignment library sources of [lh3/bwa](https://github.com/lh3/bwa) at commit
`cb950614` (0.7.16a), the commit wrapped by the gatk-bwamem-jni library that Java ESVEE uses, together
with `esvee_bwa.c`, which is the non-JNI part of gatk-bwamem-jni's `jnibwa.c`. BWA is
licensed under GPLv3; `ksw.c` and the klib headers are MIT, and so is `sse2neon.h`. A binary built from
this crate therefore falls under GPLv3, which is also the licence of hmftools ([../LICENSE](../LICENSE)).
See [vendor/bwa/README.md](vendor/bwa/README.md).
