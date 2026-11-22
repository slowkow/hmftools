# esvee-rs

A Rust port of [ESVEE](../esvee/README.md) 2.1, the structural variant caller. One binary holds the four
steps of the Java tool (prep, assembly, depth annotation, calling); each step takes the command line of
the Java application of that name and writes the same files.

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

## What is ported

| Step | Ported | Not ported |
|---|---|---|
| prep | everything: tumor/normal, tumor only, normal only, BAM and CRAM input, Illumina, SBX and Ultima, the cache BAM and `-no_cache_bam` modes, `-write_types`, SAGA location matching (`-saga_fasta`) | nothing known |
| assembly | everything that writes the TSVs and the raw VCF: tumor/normal, tumor only, germline only, Illumina, SBX and Ultima, SAGA, `-decoy_genome`, all `-write_types`, `-perf_debug` and `-perf_log_time` log lines | the assembly visualiser: `-run_assembly_vis` and `-purple_dir` are accepted and the main outputs are unchanged, but no HTML is written |
| depth | everything, BAM and CRAM (3.0 and 3.1) | htsjdk's record validation behind `-bam_validation STRICT` and `LENIENT`: the option is accepted, records are not validated |
| caller | everything: tumor/normal, tumor only, germline only, PON and artefact PON, hotspots, repeat mask, target regions, `-write_breakend_tsv`, `-manual_ref_depth` | nothing known |

What has and has not been tested for each is in REPORT.md, sections 2 and 6. In short: all four steps
give Java's output on a small Illumina tumor/normal input and across the option cases; the whole-genome
run has been made with the port but not yet compared with Java; reference genome 37 has not been run on
real data in any step.

## What "the same" means

- **BAM:** the same records in the same order with the same fields and attributes, and the same header
  apart from `@PG` lines, compared after decompression.
- **TSV:** byte-identical.
- **VCF:** byte-identical after decompression, header included.
- **Compressed bytes and index files** (`.bai`, `.tbi`) differ from Java's, since the compressor and
  the index writer differ. The indexes were checked by querying through them.
- **Thread count:** the port writes the same output at any thread count. Java does not:
  - Java's prep output depends on the thread count (the partitions share a cache of spanning reads).
    The port gives Java's one-thread result.
  - Java's assembly output depends on the thread count and also changes from run to run at one thread,
    because sets of objects are iterated in identity-hash order. The port gives the result of Java at
    one thread with every identity hash made equal (`-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`).
  - Java's depth and caller VCFs do not depend on the thread count, and the port matches them.

  To compare the port with Java, run Java with `-threads 1` and, for the assembly step, with that JVM
  option. REPORT.md, section 3, has the detail.

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
- **Faults in Java are reproduced** where they change the output, so that the output can be compared.
  They are listed in REPORT.md, section 4, and marked in the code.
- **A missing SAGA index image** stops the assembly step at start-up; Java only fails when the first
  matcher is created.
- **Log lines** have the Java text. Timings in them differ, and with `-perf_log_time` fewer items
  exceed a given threshold. The caller prints the lines about its reference data together once the
  data is loaded.
- **Memory** is lower in every step on the small input (REPORT.md, section 5). The assembly step's
  memory on a whole genome is high (33 GB) and is being worked on.

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
