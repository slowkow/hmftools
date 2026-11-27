#!/usr/bin/env python3
"""Writes a mutation contexts table of random synthetic samples (no real data): random_contexts.py
<rows> <seed> <out.tsv>. The counts are drawn from mixtures that reach all of the model's regions:
HR-proficient-like, BRCA1-like (small duplications, microhomology deletions), BRCA2-like (deletions
with microhomology, small deletions), MSI-like (repeat indels), and rows with few or no mutations of a
kind, so that the QC rules and ties are exercised too."""
import random, sys

def names():
    out = []
    for r in "CT":
        for a in "ACGT":
            if a == r:
                continue
            for b in "ACGT":
                for c in "ACGT":
                    out.append(f"{b}[{r}>{a}]{c}")
    for ctx, cap in ((".rep.len.", 5), (".mh.bimh.", 5), (".none.len.", 5)):
        for m in ("del", "ins"):
            for i in range(1, cap + 1):
                out.append(f"{m}{ctx}{i}")
    names_sv = []
    for t in ("DEL", "DUP", "INV"):
        bins = ["0e00", "1e03", "1e04", "1e05", "1e06", "1e07", "Inf"]
        for i in range(6):
            names_sv.append(f"{t}_{bins[i]}_{bins[i+1]}_bp")
    names_sv.append("TRA")
    return out + names_sv

def draw(rng, n, weights):
    tot = sum(weights)
    if n == 0 or tot == 0:
        return [0] * len(weights)
    counts = [0] * len(weights)
    cum = []
    s = 0
    for w in weights:
        s += w
        cum.append(s)
    for _ in range(n):
        u = rng.random() * tot
        lo, hi = 0, len(cum) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cum[mid] < u:
                lo = mid + 1
            else:
                hi = mid
        counts[lo] += 1
    return counts

def main():
    rows, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    rng = random.Random(seed)
    cols = names()
    with open(out, "w") as f:
        f.write("sample_id\t" + "\t".join(cols) + "\n")
        for r in range(rows):
            kind = rng.choice(["prof", "brca1", "brca2", "msi", "few", "mixed", "mixed"])
            scale = rng.choice([0, 1, 5, 30, 100, 1000, 5000, 20000])
            n_snv = int(rng.expovariate(1.0) * scale * 3)
            n_indel = int(rng.expovariate(1.0) * scale / 3)
            n_sv = int(rng.expovariate(1.0) * scale / 30)
            if kind == "few":
                n_snv, n_indel, n_sv = rng.choice([0, 1, 3]), rng.choice([0, 2, 50, 99, 100]), rng.choice([0, 1, 29, 30])
            if kind == "msi":
                n_indel = rng.randint(14000, 30000)
            snv_w = [rng.random() ** rng.choice([1, 3, 8]) for _ in range(96)]
            indel_w = [rng.random() ** rng.choice([1, 3, 8]) for _ in range(30)]
            sv_w = [rng.random() ** rng.choice([1, 3, 8]) for _ in range(19)]
            if kind == "brca1":
                for i in (7, 8, 9):  # DUP 1e03-1e06
                    sv_w[i] *= 20
                for i in range(10, 15):
                    indel_w[i] *= 10
            if kind == "brca2":
                for i in (0, 1, 2):
                    sv_w[i] *= 20
                for i in range(11, 15):
                    indel_w[i] *= 20
            if kind == "msi":
                for i in range(0, 10):
                    indel_w[i] *= 50
            counts = draw(rng, n_snv, snv_w) + draw(rng, n_indel, indel_w) + draw(rng, n_sv, sv_w)
            f.write(f"S{r + 1}\t" + "\t".join(map(str, counts)) + "\n")

main()
