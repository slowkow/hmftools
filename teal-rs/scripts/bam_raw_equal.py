#!/usr/bin/env python3
"""usage: bam_raw_equal.py <a.bam> <b.bam>: whether the uncompressed header and record bytes are identical
(every field, the bin and the attribute types included); prints the first difference"""
import gzip, sys
a = gzip.open(sys.argv[1]).read()
b = gzip.open(sys.argv[2]).read()
if a == b:
    print(f"identical bytes ({len(a)} uncompressed)")
    sys.exit(0)
i = next((k for k in range(min(len(a), len(b))) if a[k] != b[k]), min(len(a), len(b)))
print(f"differ at uncompressed byte {i} of {len(a)}/{len(b)}")
sys.exit(1)
