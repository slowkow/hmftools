# esvee-rs: status

A port of ESVEE 2.1, one module per step, each checked against the Java step's output given the same
input. The Java reference is upstream master at commit `ca04e6c39d` (2026-11-20) (the fix for prep BAM loading is needed).

Since 2026-11-22 the port corrects 61 of the Java findings by default; `-java_compat` gives Java's
behaviour, and every comparison with Java in the table below is made in that mode (REPORT.md, 4.5).

How to build and run: [README.md](README.md). What was compared, the Java findings and the
measurements in full: [REPORT.md](REPORT.md).

| Step | Status |
|---|---|
| prep | ported; prep BAMs (same records, same order) and the three TSVs identical to Java's one-thread output at any thread count, on a small tumor/normal input, an exome pair (apart from five tied record pairs that Java itself orders differently from run to run) and 42 option checks incl. public SBX and Ultima reads; on a whole genome the three TSVs identical and both BAMs the same records, with one such tied pair in each in the other order; all 21 Java unit tests ported |
| assembly | ported; raw VCF, assembly.tsv, breakend.tsv and the optional TSVs identical to Java run at one thread with constant identity hashes (`-XX:hashCode=2`), at any Rust thread count, in 52 option cases (Illumina, public SBX and Ultima, exome, germline-only, tumor-only, SAGA, decoy); BWA-MEM is the same C source as Java's JNI library, built with the crate; all 20 Java test classes ported; the assembly visualiser writes the same HTML files as Java (96 files over the 7 cases in which Java writes any, byte for byte); genome 37 untested |
| depth | ported; identical VCF to Java on a small tumor/normal input and 22 synthetic inputs, at 1 and 8 threads; not ported: htsjdk record validation behind `-bam_validation STRICT` |
| caller | ported; unfiltered, somatic and germline VCFs identical to Java on the small input and in 96 option cases; all 18 Java unit tests ported; SBX, Ultima and genome 37 not yet run on real data of those kinds |
| end to end | on the small input every output of the Rust chain at 8 threads equals the Java reference chain (one thread, `-XX:hashCode=2`): prep BAMs and TSVs, assembly.tsv, breakend.tsv, raw / ref_depth / unfiltered / somatic / germline VCFs |
| one call (EsveeApplication) | `esvee-rs <options>` without a step name runs the four steps with the options of all four, as oncoanalyser runs ESVEE; on the small input with oncoanalyser's command line (sambamba, its write types, SAGA and decoy) every output equals Java's EsveeApplication, tumor/normal and tumor only (REPORT.md, 2.5) |
| whole genome | Illumina tumor/normal pair (4.2 billion reads): assembly.tsv, breakend.tsv and all five VCFs (52,061 variants, 3,072 somatic, 35,338 germline) identical to the Java reference chain, from the Rust prep output and from Java's; prep as above; Rust 14 min at 8 threads against 2.9 h for the Java reference at one thread |

## Measurements

Prep, small input (120 M reads, two BAMs): Java 55–58 s / 40–44 GB at 8 threads and about 240 s / 5.8 GB
at one; Rust 10–22 s / 2.3–2.6 GB and 67 s / 1.4 GB. Java's prep output depends on thread count (the
partitions share a cache of spanning reads); the port gives Java's one-thread result at any thread count.

Assembly, small input: Java 6.5–18 s / 6–8 GB at 8 threads and 30 s at one; Rust 2.7–3 s / 1.5 GB and
17 s. Java's assembly output changes from run to run even at one thread (sets of objects iterated in
identity-hash order) and with thread count; the reference is Java with every identity hash made equal.

Depth, small input (2,425 variants, two BAMs): Java 10.4 s / 1.6–4.6 GB at 8 threads and 19.6 s at one;
Rust 0.9 s / 63 MB and 5.8 s / 34 MB.

Caller, small input (2,425 records): Java 11 s / 2.5–5 GB; Rust 1.1 s / 240 MB.

Whole genome (tumor/normal pair), Rust at 8 threads: prep 672 s / 5.1 GB, assembly 110-127 s / 7.2 GB,
depth 55-60 s / 0.7 GB, caller 2.6 s / 0.5 GB. Java reference at one thread: prep 8,466 s / 14.1 GB,
assembly 1,679 s / 30.4 GB, depth 413 s / 12.6 GB, caller 13.5 s / 6.5 GB. The port's assembly peak
was 33.4 GB at first and is 7.2 GB now (6.7 GB at one thread, 803 s), with the same output.
