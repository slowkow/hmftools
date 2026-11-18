#!/usr/bin/env python3
"""Compare the per-sample genotype columns of a Java SAGE VCF with the Rust -log_evidence TSV."""
import gzip, sys
from collections import Counter

java_vcf, rust_tsv = sys.argv[1], sys.argv[2]
keys = None
java = {}
with gzip.open(java_vcf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        keys = v[8].split(":")
        java.setdefault((v[0], v[1], v[3], v[4]), []).append(v[9:])
rust = {}
with open(rust_tsv) as f:
    for line in f:
        v = line.rstrip("\n").split("\t")
        rust.setdefault((v[0], v[1], v[2], v[3]), []).append(v[4:])

only_java = [k for k in java if k not in rust]
only_rust = [k for k in rust if k not in java]
field_diffs = Counter()
examples = {}
compared = 0
for k in java:
    if k not in rust:
        continue
    for jrec, rrec in zip(java[k], rust[k]):
        compared += 1
        for s, (js, rs) in enumerate(zip(jrec, rrec)):
            for name, a, b in zip(keys, js.split(":"), rs.split(":")):
                if a != b:
                    field_diffs[(s, name)] += 1
                    examples.setdefault((s, name), (k, a, b))
print(f"java records {sum(map(len, java.values()))}, rust records {sum(map(len, rust.values()))}, compared {compared}")
print(f"only in java {len(only_java)}, only in rust {len(only_rust)}")
for k in only_java[:5]:
    print("  java only", k)
for k in only_rust[:5]:
    print("  rust only", k)
for (s, name), n in sorted(field_diffs.items(), key=lambda x: -x[1]):
    k, a, b = examples[(s, name)]
    print(f"  sample {s} {name}: {n} differ, e.g. {k[0]}:{k[1]} {k[2]}>{k[3]} java={a} rust={b}")
if not field_diffs and not only_java and not only_rust:
    print("IDENTICAL")
