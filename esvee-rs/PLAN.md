# esvee-rs: status

A port of ESVEE 2.1, one module per step, each checked against the Java step's output given the same
input. The Java reference is upstream master of 2026-11-19 (the fix for prep BAM loading is needed).

| Step | Status |
|---|---|
| prep | in progress |
| assembly | in progress |
| depth | ported; identical VCF to Java on a small tumor/normal input and 22 synthetic inputs, at 1 and 8 threads; not ported: htsjdk record validation behind `-bam_validation STRICT` |
| caller | in progress |

Depth, small input (2,425 variants, two BAMs): Java 10.4 s / 1.6–4.6 GB at 8 threads and 19.6 s at one;
Rust 0.9 s / 63 MB and 5.8 s / 34 MB.
