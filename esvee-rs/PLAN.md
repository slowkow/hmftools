# esvee-rs: status

A port of ESVEE 2.1, one module per step, each checked against the Java step's output given the same
input. The Java reference is upstream master of 2026-11-19 (the fix for prep BAM loading is needed).

| Step | Status |
|---|---|
| prep | ported; prep BAMs (same records, same order) and the three TSVs identical to Java's one-thread output at any thread count, on a small tumor/normal input, an exome pair (apart from five tied record pairs that Java itself orders differently from run to run) and 42 option checks incl. public SBX and Ultima reads; all 21 Java unit tests ported; whole genome not yet compared with Java |
| assembly | in progress |
| depth | ported; identical VCF to Java on a small tumor/normal input and 22 synthetic inputs, at 1 and 8 threads; not ported: htsjdk record validation behind `-bam_validation STRICT` |
| caller | ported; unfiltered, somatic and germline VCFs identical to Java on the small input and in 96 option cases; all 18 Java unit tests ported; SBX, Ultima and genome 37 not yet run on real data of those kinds |

Depth, small input (2,425 variants, two BAMs): Java 10.4 s / 1.6–4.6 GB at 8 threads and 19.6 s at one;
Rust 0.9 s / 63 MB and 5.8 s / 34 MB.

Caller, small input (2,425 records): Java 11 s / 2.5–5 GB; Rust 1.1 s / 240 MB.

Prep, small input (120 M reads, two BAMs): Java 55–58 s / 40–44 GB at 8 threads and about 240 s / 5.8 GB
at one; Rust 10–22 s / 2.3–2.6 GB and 67 s / 1.4 GB. Java's prep output depends on thread count (the
partitions share a cache of spanning reads); the port gives Java's one-thread result at any thread count.
