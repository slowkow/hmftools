#!/usr/bin/env python3
"""Add simulated tumor-in-normal contamination to a SAGE VCF, to exercise TINC detection and recovery.

For each variant, Binomial(reference depth, level * tumor AF) alt fragments are moved from the reference
sample's ref support to its alt support, and germline soft filters are added to the variants that gained
alt support. The filters are assigned at random rather than recomputed, so the output is only useful for
comparing two TINC implementations on identical input, not as a realistic call set.

usage: simulate_tinc.py <in.vcf.gz> <out.vcf.gz> <level> <reference_id> <tumor_id> [seed]
"""
import gzip
import random
import subprocess
import sys

GERMLINE = [("maxGermlineVAF", 0.8), ("maxGermlineRelQual", 0.6), ("maxGermlineAltSupport", 0.3)]


def binomial(rng, n, p):
    return sum(1 for _ in range(n) if rng.random() < p)


def main():
    src, dst, level, ref_id, tumor_id = sys.argv[1:6]
    level = float(level)
    rng = random.Random(int(sys.argv[6]) if len(sys.argv) > 6 else 1)
    out = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(dst, "wb"))
    ref_col = tumor_col = None
    changed = 0
    with gzip.open(src, "rt") as f:
        for line in f:
            if line.startswith("##"):
                out.stdin.write(line.encode())
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                ref_col, tumor_col = fields.index(ref_id), fields.index(tumor_id)
                out.stdin.write(line.encode())
                continue
            keys = fields[8].split(":")
            ref = dict(zip(keys, fields[ref_col].split(":")))
            tumor = dict(zip(keys, fields[tumor_col].split(":")))
            tumor_dp = int(tumor["DP"])
            tumor_af = int(tumor["AD"].split(",")[1]) / tumor_dp if tumor_dp > 0 else 0.0
            added = binomial(rng, int(ref["DP"]), min(level * tumor_af, 1.0))
            if added > 0:
                changed += 1
                ad = [int(x) for x in ref["AD"].split(",")]
                ad = [max(ad[0] - added, 0), ad[1] + added]
                ref["AD"] = ",".join(map(str, ad))
                counts = [int(x) for x in ref["RC_CNT"].split(",")]
                counts[0] += added
                counts[4] = max(counts[4] - added, 0)
                ref["RC_CNT"] = ",".join(map(str, counts))
                quals = [int(x) for x in ref["RC_QUAL"].split(",")]
                quals[0] += added * rng.randint(20, 37)
                ref["RC_QUAL"] = ",".join(map(str, quals))
                abq = ref["ARCBQ"].split(",")
                if abq[1] == "0":
                    abq[1] = str(rng.randint(25, 37))
                ref["ARCBQ"] = ",".join(abq)
                fields[ref_col] = ":".join(ref[k] for k in keys)
                filters = set(x for x in fields[6].split(";") if x not in ("PASS", "."))
                picked = [name for name, p in GERMLINE if rng.random() < p] or [GERMLINE[0][0]]
                filters.update(picked)
                fields[6] = ";".join(sorted(filters))
            out.stdin.write(("\t".join(fields) + "\n").encode())
    out.stdin.close()
    if out.wait() != 0:
        sys.exit("bgzip failed")
    subprocess.check_call(["tabix", "-f", "-p", "vcf", dst])
    print(f"variants given extra reference alt support: {changed}", file=sys.stderr)


if __name__ == "__main__":
    main()
