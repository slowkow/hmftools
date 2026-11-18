#!/usr/bin/env python3
"""Compare two SAGE VCFs record by record.

Java SAGE is not reproducible across thread counts in two INFO fields: LPS ids come from a counter
shared between threads, and NEAR_HOTSPOT is set after the record is queued for writing. This tool
reports differences in those two fields separately from everything else.
"""
import gzip, sys
from collections import Counter

def read(path):
    recs = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                recs.append(line.rstrip("\n").split("\t"))
    return recs

def split_info(info):
    d = {}
    for item in info.split(";"):
        k, _, v = item.partition("=")
        d[k] = v
    return d

a, b = read(sys.argv[1]), read(sys.argv[2])
key = lambda r: (r[0], int(r[1]), r[3], r[4])
ca, cb = Counter(map(key, a)), Counter(map(key, b))
only_a = sum((ca - cb).values())
only_b = sum((cb - ca).values())
print(f"records: {len(a)} vs {len(b)}; only in first {only_a}, only in second {only_b}")
for k in list((ca - cb).keys())[:5]:
    print("  only in first:", k)
for k in list((cb - ca).keys())[:5]:
    print("  only in second:", k)

diffs = Counter()
examples = {}
ia = ib = 0
compared = identical = 0
# both files are in the same order; walk them together, skipping records unique to one side
while ia < len(a) and ib < len(b):
    ka, kb = key(a[ia]), key(b[ib])
    if ka != kb:
        if cb[ka] == 0 or (ka[1] < kb[1] and ka[0] == kb[0]):
            ia += 1
        else:
            ib += 1
        continue
    ra, rb = a[ia], b[ib]
    ia += 1
    ib += 1
    compared += 1
    if ra == rb:
        identical += 1
        continue
    def note(name):
        diffs[name] += 1
        examples.setdefault(name, f"{ra[0]}:{ra[1]} {ra[3]}>{ra[4]}")
    if ra[5] != rb[5]:
        note("QUAL")
    if ra[6] != rb[6]:
        note("FILTER")
    fa, fb = split_info(ra[7]), split_info(rb[7])
    for k in sorted(set(fa) | set(fb)):
        if fa.get(k) != fb.get(k):
            if k == "LPS" and k in fa and k in fb and len(fa[k].split(",")) == len(fb[k].split(",")):
                note("INFO/LPS (ids only)")
            else:
                note("INFO/" + k)
    for s, (ga, gb) in enumerate(zip(ra[9:], rb[9:])):
        if ga != gb:
            for name, x, y in zip(ra[8].split(":"), ga.split(":"), gb.split(":")):
                if x != y:
                    note(f"sample{s}/{name}")
print(f"compared {compared}, identical lines {identical}")
for name, n in diffs.most_common():
    print(f"  {name}: {n}   e.g. {examples[name]}")
