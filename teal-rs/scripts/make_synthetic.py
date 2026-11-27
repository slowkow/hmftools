#!/usr/bin/env python3
"""Synthetic tumor and normal BAMs for comparing Java TEAL with teal-rs (no real data).

usage: make_synthetic.py <out dir> [seed] [unknown]   (unknown: some SA tags name a contig not in the header)

Writes <out>/{tumor,normal}.sam (coordinate order is left to samtools sort) plus the pipeline files the
TealPipelineApp needs (cobalt/, purple/, bamtools/). Contigs chr1-chr3 (5 Mb), chrX (2 Mb) and one alt
contig. Fragments: background pairs, telomeric pairs with mates near and far, split reads with telomeric
clips at shared (germline) and tumor-only break end sites (clusters within a few bases), supplementary
alignments with SA tags (optionally some to a contig not in the header), unplaced unmapped pairs in both mate orders, unmapped
placed mates, secondary and consensus (CR) records, duplicates, reads with N, and single-end reads.
"""
import os
import random
import sys

UNKNOWN_CONTIG = False
CONTIGS = [("chr1", 5_000_000), ("chr2", 5_000_000), ("chr3", 5_000_000), ("chrX", 2_000_000), ("chr1_KI270706v1_random", 200_000)]
G = "TTAGGG"
C = "CCCTAA"
READ_LEN = 150


def rnd_seq(r, n):
    return "".join(r.choice("ACGT") for _ in range(n))


def telo(r, n, g_rich=True, errors=0.02, offset=None):
    hexamer = G if g_rich else C
    off = r.randrange(6) if offset is None else offset
    s = (hexamer * (n // 6 + 3))[off:off + n]
    out = []
    for ch in s:
        if r.random() < errors:
            out.append(r.choice("ACGT"))
        else:
            out.append(ch)
    return "".join(out)


def revcomp(s):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(comp[c] for c in reversed(s))


class Sample:
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.n = 0

    def qname(self, kind):
        self.n += 1
        return f"SYN:{kind}:{self.n:06d}"

    def add(self, qname, flag, chrom, pos, mapq, cigar, rnext, pnext, tlen, seq, tags=()):
        qual = "".join(chr(33 + 30) for _ in seq) if seq != "*" else "*"
        self.lines.append("\t".join([qname, str(flag), chrom, str(pos), str(mapq), cigar, rnext, str(pnext), str(tlen), seq, qual] + list(tags)))

    def pair(self, qname, c1, p1, s1, cig1, c2, p2, s2, cig2, rev1=False, rev2=True, extra1=(), extra2=(), dup=False, mapq=60, flags1=0, flags2=0):
        """a pair; c = None for an unmapped read (placed at its mate when the mate is mapped)"""
        f1 = 1 | 64 | (16 if rev1 else 0) | (32 if rev2 else 0) | (1024 if dup else 0) | flags1
        f2 = 1 | 128 | (16 if rev2 else 0) | (32 if rev1 else 0) | (1024 if dup else 0) | flags2
        if c1 is None:
            f1 |= 4
            f2 |= 8
        if c2 is None:
            f2 |= 4
            f1 |= 8
        # unmapped reads take the mate's position, or none
        r1c, r1p = (c1, p1) if c1 is not None else ((c2, p2) if c2 is not None else ("*", 0))
        r2c, r2p = (c2, p2) if c2 is not None else ((c1, p1) if c1 is not None else ("*", 0))
        mn1 = "=" if r2c == r1c and r1c != "*" else r2c
        mn2 = "=" if r1c == r2c and r2c != "*" else r1c
        tlen = 0
        if c1 is not None and c1 == c2:
            tlen = (p2 + READ_LEN) - p1 if p2 >= p1 else -((p1 + READ_LEN) - p2)
        self.add(qname, f1, r1c, r1p, mapq if c1 else 0, cig1 if c1 else "*", mn1, r2p, tlen, s1, extra1)
        self.add(qname, f2, r2c, r2p, mapq if c2 else 0, cig2 if c2 else "*", mn2, r1p, -tlen, s2, extra2)


def contig(r, exclude_alt=False):
    names = CONTIGS[:-1] if exclude_alt else CONTIGS
    weights = [l for _, l in names]
    return r.choices(names, weights=weights)[0]


def pos_in(r, c, margin=1000):
    return r.randrange(margin, c[1] - margin)


def build(sample, r, sites, tumor):
    # background pairs
    for _ in range(6000):
        c = contig(r)
        p = pos_in(r, c)
        sample.pair(sample.qname("bg"), c[0], p, rnd_seq(r, READ_LEN), "150M", c[0], p + r.randrange(100, 400), rnd_seq(r, READ_LEN), "150M")
    # telomeric pairs: the telomeric read before or after its mate, mates near or on another contig
    for i in range(400):
        c = contig(r)
        p = pos_in(r, c)
        g_rich = r.random() < 0.5
        tseq = telo(r, READ_LEN, g_rich)
        if r.random() < 0.6:
            c2, p2 = c, p + r.choice([-1, 1]) * r.randrange(200, 600)
        else:
            c2 = contig(r)
            p2 = pos_in(r, c2)
        mseq = rnd_seq(r, READ_LEN) if r.random() < 0.7 else telo(r, READ_LEN, not g_rich)
        q = sample.qname("tp")
        if r.random() < 0.5:
            sample.pair(q, c[0], p, tseq, "150M", c2[0], p2, mseq, "150M")
        else:
            sample.pair(q, c2[0], p2, mseq, "150M", c[0], p, tseq, "150M", rev1=True, rev2=False)
    # split reads at break end sites: telomeric clip + anchored part, mate telomeric (mapped near or
    # unmapped placed), some with supplementary alignments
    for (sc, spos, kind, count, germline) in sites:
        if not tumor and not germline:
            continue
        n = count if tumor else max(1, count // 3)
        for k in range(n):
            pos = spos + r.choice([0, 0, 0, 1, 2, -1, 3])
            clip = r.randrange(20, 80)
            anchor = READ_LEN - clip
            g_rich = kind in ("RIGHT_G", "LEFT_G")
            q = sample.qname("sr")
            dup = r.random() < 0.1
            if kind.startswith("LEFT"):
                seq = telo(r, clip, g_rich) + rnd_seq(r, anchor)
                cig = f"{clip}S{anchor}M"
                start = pos
                rev = True
            else:
                seq = rnd_seq(r, anchor) + telo(r, clip, g_rich)
                cig = f"{anchor}M{clip}S"
                start = pos - anchor + 1
                rev = False
            mate = telo(r, READ_LEN, g_rich if r.random() < 0.7 else not g_rich)
            tags1 = []
            supp = None
            if r.random() < 0.3:
                # the clip aligned elsewhere (as supplementary), sometimes on a contig not in the header
                if UNKNOWN_CONTIG and r.random() < 0.2:
                    sa_c, sa_p = "chrUnknown", 1234
                else:
                    cc = contig(r, exclude_alt=True)
                    sa_c, sa_p = cc[0], pos_in(r, cc)
                sa_cig = f"{anchor}S{clip}M" if kind.startswith("LEFT") else f"{anchor}S{clip}M"
                tags1.append(f"SA:Z:{sa_c},{sa_p},+,{sa_cig},20,1;")
                if sa_c != "chrUnknown":
                    supp = (sa_c, sa_p, sa_cig, f"SA:Z:{sc},{start},{'-' if rev else '+'},{cig},60,0;")
            if r.random() < 0.25:
                # mate unmapped, placed at the split read
                sample.pair(q, sc, start, seq, cig, None, 0, mate, "*", rev1=rev, rev2=not rev, extra1=tags1, dup=dup)
            else:
                mpos = start + (r.randrange(150, 400) if not rev else -r.randrange(150, 400))
                sample.pair(q, sc, start, seq, cig, sc, max(1, mpos), mate, "150M", rev1=rev, rev2=not rev, extra1=tags1, dup=dup)
            if supp:
                f = 1 | 64 | 2048 | (32 if not rev else 0)
                sample.add(q, f, supp[0], supp[1], 20, supp[2], sc, start, 0, seq, [supp[3]])
    # unplaced unmapped pairs, the telomeric read first or second, one or both telomeric
    for _ in range(150):
        q = sample.qname("uu")
        a = telo(r, READ_LEN, r.random() < 0.5)
        b = rnd_seq(r, READ_LEN) if r.random() < 0.5 else telo(r, READ_LEN, r.random() < 0.5)
        if r.random() < 0.5:
            a, b = b, a
        sample.add(q, 1 | 4 | 8 | 64, "*", 0, 0, "*", "*", 0, 0, a)
        sample.add(q, 1 | 4 | 8 | 128, "*", 0, 0, "*", "*", 0, 0, b)
    # an unplaced unmapped telomeric read with an SA attribute, and one whose mate is missing
    q = sample.qname("uusa")
    sample.add(q, 1 | 4 | 8 | 64, "*", 0, 0, "*", "*", 0, 0, telo(r, READ_LEN, True), ["SA:Z:chr2,3000000,+,50S100M,0,1;"])
    sample.add(q, 1 | 4 | 8 | 128, "*", 0, 0, "*", "*", 0, 0, rnd_seq(r, READ_LEN))
    sample.add(sample.qname("uulone"), 1 | 4 | 8 | 64, "*", 0, 0, "*", "*", 0, 0, telo(r, READ_LEN, False))
    # secondary alignments: telomeric, of a complete pair (they open a new group) and of a background pair
    for _ in range(20):
        c = contig(r, exclude_alt=True)
        p = pos_in(r, c)
        q = sample.qname("sec")
        tseq = telo(r, READ_LEN, True)
        sample.pair(q, c[0], p, tseq, "150M", c[0], p + 300, rnd_seq(r, READ_LEN), "150M")
        c2 = contig(r, exclude_alt=True)
        sample.add(q, 1 | 64 | 256, c2[0], pos_in(r, c2), 0, "150M", c[0], p + 300, 0, tseq)
    # consensus reads (CR attribute): telomeric, skipped; their mates are not
    for _ in range(20):
        c = contig(r, exclude_alt=True)
        p = pos_in(r, c)
        q = sample.qname("cr")
        sample.pair(q, c[0], p, telo(r, READ_LEN, False), "150M", c[0], p + 250, rnd_seq(r, READ_LEN), "150M", extra1=["CR:Z:3;3;0"])
    # reads with N next to and inside telomeric repeats
    for _ in range(60):
        c = contig(r, exclude_alt=True)
        p = pos_in(r, c)
        s = list(telo(r, READ_LEN, r.random() < 0.5, errors=0.0))
        for _ in range(r.randrange(1, 6)):
            s[r.randrange(READ_LEN)] = "N"
        q = sample.qname("nn")
        sample.pair(q, c[0], p, "".join(s), "150M", c[0], p + 200, telo(r, READ_LEN, r.random() < 0.5), "150M", rev1=r.random() < 0.5, rev2=r.random() < 0.5)
    # single-end reads (no paired flag), telomeric, mapped and unmapped, some with SA
    for _ in range(60):
        c = contig(r, exclude_alt=True)
        p = pos_in(r, c)
        q = sample.qname("se")
        s = rnd_seq(r, 40) + telo(r, 110, r.random() < 0.5)
        tags = []
        if r.random() < 0.3:
            c2 = contig(r, exclude_alt=True)
            p2 = pos_in(r, c2)
            tags = [f"SA:Z:{c2[0]},{p2},+,40S110M,10,2;"]
            sample.add(q, 2048, c2[0], p2, 10, "40S110M", "*", 0, 0, s, [f"SA:Z:{c[0]},{p},+,40M110S,60,0;"])
        sample.add(q, 16 if r.random() < 0.5 else 0, c[0], p, 60, "40M110S", "*", 0, 0, s, tags)
    sample.add(sample.qname("seu"), 4, "*", 0, 0, "*", "*", 0, 0, telo(r, READ_LEN, True))
    # a pair ending at base 9,000 of chr1, at the start of the first GRCh37 blacklisted region (BED start
    # 9000): excluded when the BED start is read as 1-based (TE7)
    sample.pair(sample.qname("bl"), "chr1", 8851, telo(r, READ_LEN, True), "150M", "chr1", 9100, rnd_seq(r, READ_LEN), "150M")
    # split reads starting at base 9,000 of chr1 (a C-rich left clip): a break end candidate only when the BED
    # start is read as 0-based (TE7 corrected), since Java excludes reads lying within 9000-11000
    for _ in range(4 if tumor else 1):
        q = sample.qname("blsr")
        sample.pair(q, "chr1", 9000, telo(r, 40, False) + rnd_seq(r, 110), "40S110M", "chr1", 8700, telo(r, READ_LEN, False), "150M", rev1=True, rev2=False)
    # a telomeric pair in the poly-G region of LINC00486 is no candidate there (v38 chr2:32,916,190) - outside
    # these contigs; a pair on the alt contig instead
    c = CONTIGS[-1]
    sample.pair(sample.qname("alt"), c[0], 5000, telo(r, READ_LEN, True), "150M", c[0], 5300, rnd_seq(r, READ_LEN), "150M")


def write_inputs(out):
    os.makedirs(f"{out}/cobalt", exist_ok=True)
    os.makedirs(f"{out}/purple", exist_ok=True)
    os.makedirs(f"{out}/bamtools", exist_ok=True)
    for s, mean, gc50 in (("TUMOR", 61.25, 47.5), ("NORMAL", 33.1, 31.0)):
        with open(f"{out}/cobalt/{s}.cobalt.gc.median.tsv", "w") as f:
            f.write(f"#sampleMean\tsampleMedian\n{mean}\t{mean - 1}\n#gcBucket\tmedian\n")
            for b in range(20, 70):
                f.write(f"{b}\t{gc50 if b == 50 else mean + b / 10:.2f}\n")
        with open(f"{out}/bamtools/{s}.bam_metric.summary.tsv", "w") as f:
            cols = "TotalRegionBases\tTotalReads\tDuplicateReads\tDualStrandReads\tOffTargetReads\tMeanCoverage\tStdDevCoverage\tMedianCoverage\tMadCoverage\tLowMapQualPercent\tDuplicatePercent\tUnmappedPercent\tLowBaseQualPercent\tOverlappingReadPercent\tCappedCoverage\tDepthCoverage_1\tDepthCoverage_10"
            dup = "0.13757" if s == "NORMAL" else "0.12275"
            f.write(cols + "\n" + f"1000000\t20000\t2000\t0\t0\t30.000\t10.000\t30\t5\t0.04000\t{dup}\t0.00100\t0.00030\t0.22000\t0.00010\t0.97\t0.90\n")
    with open(f"{out}/purple/TUMOR.purple.purity.tsv", "w") as f:
        f.write("purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\trunMode\ttargeted\n")
        f.write("0.4200\t0.7000\t0.7000\t0.1000\t2.6500\tFEMALE\tNORMAL\t0.1000\t0.4000\t0.4500\t2.5000\t2.8000\t0.0500\t0.1500\t0.0000\tfalse\t0.0300\tMSS\t10\tLOW\t1.0000\tLOW\t50\tTUMOR_GERMLINE\tfalse\n")
    with open(f"{out}/purple/TUMOR.purple.qc", "w") as f:
        f.write("QCStatus\tPASS\nMethod\tNORMAL\nCopyNumberSegments\t500\nUnsupportedCopyNumberSegments\t10\nPurity\t0.4200\nAmberGender\tFEMALE\nCobaltGender\tFEMALE\nDeletedGenes\t0\nContamination\t0.0\nGermlineAberrations\tNONE\nAmberMeanDepth\t40\nLohPercent\t0.1\nTincLevel\t0.0\n")


def main():
    out = sys.argv[1]
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    global UNKNOWN_CONTIG
    UNKNOWN_CONTIG = len(sys.argv) > 3 and sys.argv[3] == "unknown"
    r = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    # break end sites: (contig, position, type, split reads in the tumor, also germline)
    sites = []
    for i in range(40):
        c = contig(r, exclude_alt=True)
        p = pos_in(r, c, 5000)
        kind = r.choice(["RIGHT_G", "LEFT_C", "RIGHT_C", "LEFT_G"])
        sites.append((c[0], p, kind, r.randrange(1, 9), r.random() < 0.4))
        if r.random() < 0.3:
            # a second site close by (duplicate marking, merging)
            sites.append((c[0], p + r.randrange(3, 70), kind, r.randrange(1, 6), r.random() < 0.5))
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{n}\tLN:{l}\n" for n, l in CONTIGS) + "@RG\tID:syn\tSM:syn\n"
    for name, tumor in (("tumor", True), ("normal", False)):
        s = Sample(name)
        build(s, r, sites, tumor)
        with open(f"{out}/{name}.sam", "w") as f:
            f.write(header)
            for line in s.lines:
                f.write(line + "\n")
    write_inputs(out)


if __name__ == "__main__":
    main()
