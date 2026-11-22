# Vendored BWA-MEM

The ESVEE assembly step aligns assembled sequences with BWA-MEM. The Java tool does this through the
gatk-bwamem-jni library, which wraps lh3/bwa at commit cb950614ce7217788780b9a8d445c64cd4d8f62e
(0.7.16a-r1181) and reads the genome from a BWA index image (`<fasta>.img`). To produce the same alignment
records, this directory holds:

- the sources of that bwa commit that make up `libbwa.a` (the alignment library, without the index builder);
- `sse2neon.h` (DLTcollab/sse2neon v1.7.0, MIT), included by `ksw.c` in place of `<emmintrin.h>` on ARM so the
  SSE2 Smith-Waterman code builds there unchanged (the results do not depend on the SIMD implementation);
- `esvee_bwa.c`: the non-JNI part of gatk-bwamem-jni's `jnibwa.c` (open an index image, align a batch of
  sequences into the same packed result buffer), plus a function that fills in the alignment options.

BWA is licensed under GPLv3 (`COPYING`); `ksw.c` and the klib headers are MIT.
