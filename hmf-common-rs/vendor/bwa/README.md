# Vendored BWA-MEM

ESVEE (assembled sequences) and CIDER (VDJ sequences) align with BWA-MEM. The Java tools do this through the
gatk-bwamem-jni library, which wraps lh3/bwa at commit cb950614ce7217788780b9a8d445c64cd4d8f62e
(0.7.16a-r1181) and reads the genome from a BWA index image (`<fasta>.img`). To produce the same alignment
records, this directory holds:

- the sources of that bwa commit that make up `libbwa.a` (the alignment library, without the index builder);
- `sse2neon.h` (DLTcollab/sse2neon v1.7.0, MIT), included by `ksw.c` in place of `<emmintrin.h>` on ARM so the
  SSE2 Smith-Waterman code builds there unchanged (the results do not depend on the SIMD implementation);
- `hmf_bwa.c`: the non-JNI part of gatk-bwamem-jni's `jnibwa.c` (open an index image, from a file or from
  memory, align a batch of sequences into the same packed result buffer), and access to the `mem_opt_t` of
  `mem_opt_init()`, which the Rust side changes field by field as the Java aligners do. Copied from esvee-rs's
  `esvee_bwa.c` (same code, `hmf_bwa_` names).

The library is built only with the crate's `bwa` feature (see `build.rs` and `src/bwa.rs`).

BWA is licensed under GPLv3 (`COPYING`); `ksw.c` and the klib headers are MIT.
