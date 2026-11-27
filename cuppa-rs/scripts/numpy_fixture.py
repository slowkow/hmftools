"""The numbers numpy / pandas give for the operations cuppa-rs reproduces, as a test fixture.

python -I numpy_fixture.py ../testdata/numpy_fixture.tsv   (in pycuppa's environment, OPENBLAS_NUM_THREADS=1)

Lines: kind, then tab-separated fields; floats as hex (float.hex), lists comma-separated. Random inputs (seed 11).
"""
import sys

import numpy as np
import pandas as pd

rng = np.random.default_rng(11)
out = open(sys.argv[1], "w")
H = lambda v: float(v).hex()
L = lambda a: ",".join(H(v) for v in a)


def rnd(n):
    return rng.standard_normal(n) * 10.0 ** rng.integers(-3, 4, size=n)


# np.sum of a contiguous array (pairwise); the row sums of cohorts depend on the frames' layout and are
# checked by the cohort regression cases
for n in list(range(0, 20)) + [40, 96, 127, 128, 129, 200, 300, 1000, 6188]:
    a = rnd(n)
    out.write(f"sum\t{L(a)}\t{H(np.sum(a))}\n")

# dgemv_t: x @ W.T (W row-major), np.dot(x, P) (P column-major); dgemv_n: x @ W32.T (float32 W); dgemm: 2 rows
for k in [1, 3, 4, 5, 31, 32, 33, 40, 96, 120, 257, 336, 1000]:
    x = rnd(k).reshape(1, k)
    W = rnd(3 * k).reshape(3, k)
    r = (x @ W.T)[0]
    out.write(f"dot\t{L(W[1])}\t{L(x[0])}\t{H(r[1])}\n")
    W32 = W.astype(np.float32)
    r32 = (x @ W32.T)[0]
    out.write(f"dot_n\t{L(W32[2].astype(np.float64))}\t{L(x[0])}\t{H(r32[2])}\n")
    X2 = rnd(2 * k).reshape(2, k)
    g = X2 @ W.T
    out.write(f"gemm\t{L(X2[1])}\t{L(W[0])}\t{H(g[1, 0])}\n")

# np.interp, with repeated knots and values outside
for t in range(60):
    n = int(rng.integers(1, 30))
    xp = np.sort(np.round(rng.random(n), 2))
    fp = np.sort(rng.random(n))
    x = float(rng.choice([rng.random() * 1.2 - 0.1, xp[int(rng.integers(0, n))]]))
    out.write(f"interp\t{H(x)}\t{L(xp)}\t{L(fp)}\t{H(np.interp(x, xp, fp))}\n")

# argsort (default kind) with ties and NaN
for t in range(40):
    n = int(rng.integers(1, 60))
    a = rng.integers(0, 5, size=n).astype(np.float64)
    if n > 3:
        a[rng.integers(0, n)] = np.nan
    out.write(f"argsort\t{L(a)}\t{','.join(str(i) for i in np.argsort(a))}\n")

# repr (and str of a float64, as DataFrame.to_csv and astype(str) write it)
for t in range(400):
    v = float(rng.standard_normal() * 10.0 ** rng.integers(-320, 308))
    assert str(np.float64(v)) == repr(v)
    out.write(f"repr\t{H(v)}\t{repr(v)}\n")
for v in [0.0, -0.0, 1.0, 0.1, 1e16, 9999999999999998.0, 1e-4, 9.999e-5, 123456.789, 5e-324, 1.7976931348623157e308]:
    out.write(f"repr\t{H(v)}\t{repr(v)}\n")

# np.round(x, 1) and CuppaPredSummaryBuilder._signif(x, 2)
for t in range(200):
    v = float(rng.standard_normal() * 10.0 ** rng.integers(-3, 5))
    m = 10 ** (2 - 1 - np.floor(np.log10(np.abs(v))))
    out.write(f"round1\t{H(v)}\t{H(np.round(v, 1))}\n")
    out.write(f"signif2\t{H(v)}\t{H(np.round(v * m) / m)}\n")
out.close()
