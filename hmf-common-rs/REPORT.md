# hmf-common-rs: record of the ports of shared Java code

## Segmentation (2026-11-26)

hmf-common `segmentation` (PerArmSegmenter, ChrArmLocator with the GRCh37 and GRCh38 centromeres,
AbsoluteSegments, WindowSegments, SegmentsFile), `segmentation.copynumber` (Segmenter,
PiecewiseConstantFit, Gamma, Runmed, WindowedMedian, Stats: a Java port of the R copynumber
package's pcf) and `utils.pcf` (the .pcf file names, readers and the older writer). AMBER segments
its BAFs with it (absolute values, penalty factor 100), COBALT its ratios (windowed, after a float log2
taken with commons-math3 `FastMath.log`), and PURPLE reads the files.

### How Java's numbers are reproduced

- The segment costs are formed with Java's cumulative sums and operation order; ties go to the first
  start, as Java's strict `<`. The inner loop runs in eight lanes, each keeping its own first minimum, and
  merging them by cost and start gives the sequential answer.
- The running median keeps Java's two `PriorityQueue<Double>`s as sorted vectors in `Double.compare`
  order, with `remove(Object)` by `Double.equals`, so it gives Java's result with -0.0, NaN and
  infinities too.
- `DecimalFormat("#.####")` rounds the digits of the JDK's FloatingDecimal conversion rather than the
  exact binary value, so 0.00005 prints as "0" where exact rounding (Rust's `{:.4}`) gives "0.0001".
  `fmt::decimal_format` ports that conversion and DigitList's HALF_EVEN rounding.
- `FastMath.log` is ported operation by operation with its tables: `FastMath.log(2, 3)` is
  1.5849625007211563 where `Math.log(3) / Math.log(2)` is 1.584962500721156.

### Checks, in Java mode

| Check | Result |
|---|---|
| Java unit tests of these classes, ported (with the R comparisons of ExamplesTest) | 55 pass; one large-arm test runs with `--ignored` (12 s) |
| 600 synthetic cases through Java's own BAFSegmenter, a PerArmSegmenter subclass and COBALT's CobaltRatioSegmenter: random, piecewise constant, tiny and empty arms, ties, NaN, -1, -0.0, Infinity, unsorted input, both genome versions, gammas from 0 to 1e6 | 574 files byte-identical, 26 the same Java error, 0 different |
| 1 and 2 threads | identical |
| DecimalFormat digits and flags, FastMath.log, float log2 | 3,000,016 values identical |
| Running median on 20,000 series with signed zeros, NaN and infinities | identical |
| AMBER on a small tumor/normal pair: BAF segments against Java AMBER's own .pcf | byte-identical, in Java and in corrected mode |
| COBALT, whole genome: 3,088,257 windows, tumor and reference ratios, against the Java segmenter | identical |
| The .pcf readers against Java on 602 files | identical |

### Findings

Corrected by default; `-java_compat` gives Java's behaviour. A port that calls this code accepts these
ids in `-java_compat_except` (`hmf_common::compat::common_id`).

| Id | Java | Correction | Effect |
|---|---|---|---|
| SG1 | When no value is left after the `>= 0` filter, `Gamma` throws "Input array must not be empty". AMBER guards only an empty BAF list, COBALT nothing. | A header-only file. | COBALT fails on a run with no valid ratio. AMBER cannot reach it. |
| SG2 | Infinity passes the `>= 0` filter. Every cost after it is NaN, so Java either loops until the heap is exhausted (OutOfMemoryError at -Xmx256m) or writes a segment with mean "∞", which Java cannot read back. | Non-finite values are dropped like negative ones. In Java mode the never-ending case is an error. | Robustness only: neither AMBER nor COBALT produces infinite values. |
| SG3 | The running-median heaps decide with `<=` (-0.0 equals 0.0) but remove by `Double.equals` (it does not). With both zeros in one window a stale value stays: `runmed({0, -0, -0, 1, 1}, 3)` gives median 0 for the window (-0, 1, 1). | Zeros are folded to +0.0 inside the heaps. | None in AMBER or COBALT, which never pass -0.0. |

### Timings (one machine, Apple M-series)

| | Java | Rust |
|---|---|---|
| One arm of 108,041 points, 1 thread | 13.0 s, 164 MB | 3.6 s, 17 MB |
| One arm of 134,836 points, 1 thread | 16.3 s, 198 MB | 5.7 s, 19 MB |
| COBALT whole genome, tumor ratios, 2 threads | 119 s, 2.25 GB | 34.7 s, 429 MB |
| COBALT whole genome, reference ratios, 2 threads | 76 s, 2.26 GB | 34.8 s, 386 MB |

The fit is still quadratic in the points of an arm, as Java's. An exact pruned search is possible (the
costs from the stored sums are superadditive) but saves nothing on an arm without change points, such
as a flat COBALT reference.

### Remarks, not findings

- The one penalty comes from a running median over all arms joined end to end, so its windows span arm
  ends; the Java comment calls the order arbitrary.
- The fork's segmentation differs from COBALT 3.0's (3.0 used a penalty per arm for whole genomes), so
  COBALT 3.0 output is no reference for this code.
