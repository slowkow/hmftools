# hmf-common-rs

Code shared by the Rust ports of hmftools (`sage-rs`, `redux-rs`, `esvee-rs`, `amber-rs`, `cobalt-rs`,
`bam-tools-rs`, and the ports that follow).
It plays the part of the Java `hmf-common` module, plus the Java library behaviour the ports reproduce
so that their output equals the Java tools'.

| Module | What it is |
|---|---|
| `args` | command lines in the style of the Java ConfigBuilder (`-name value`, `-flag`) |
| `compat` | the switch between a port's corrected behaviour and Java's (`-java_compat`, `-java_compat_except`, `<TOOL>_JAVA_COMPAT`); each port keeps its own finding ids |
| `faidx`, `refgenome` | indexed FASTA access, 1-based inclusive positions as the Java RefGenomeInterface |
| `fmt` | numbers as Java's `String.format`, `Double.toString` and `DecimalFormat` write them (`decimal_format` ports the JDK's FloatingDecimal digits and DigitList rounding) |
| `log` | log lines in the layout of the Java tools; the `info!`, `warn!`, `error!` and `debug!` macros |
| `regions` | chromosome ranks and names, regions, `-specific_regions` / `-specific_chr`; RefGenomeVersion, BaseRegion, region and BED file loading with Java's merging, PartitionUtils |
| `java_hash`, `java_id_hash` | java.util.HashMap and HashSet iteration order (with `retain` for retainAll / removeIf, which keep the table size), for hash codes and for identity hashes under `-XX:hashCode=2` |
| `java_sort`, `java_math` | TimSort as `List.sort`, and Java numeric behaviour such as `DoubleStream.sum`, `Math.max`, `Math.round`, `Double.compare` |
| `java_fastmath` | commons-math3 `FastMath.log` and `FastMath.exp`, which differ from `Math.log` and `Math.exp` in the last bit for some inputs (COBALT takes logs with it; AMBER's Poisson and binomial distributions take both) |
| `segmentation` | hmf-common `segmentation`, `segmentation.copynumber` and `utils.pcf`: the per-arm piecewise constant fit (a port of the R copynumber package's pcf), its penalty, the segments file (AMBER absolute, COBALT windowed) and the .pcf readers (PURPLE) |
| `ensembl` | the Ensembl data cache: genes, transcripts and exons as EnsemblDataCache loads them (canonical only, restricted gene ids, synonyms), the name and id maps in Java HashMap order |
| `driver_panel` | the driver gene panel file (read, write) and DriverGeneRegions (the coding exon regions of the panel's canonical transcripts) |
| `metrics` | the BamTools metrics files: summary, coverage, fragment lengths, flag counts (BamFlagStats), gene and exon depth (GeneDepthFile); writers and readers |
| `java_text` | Java's Integer/Long/Double/Boolean parsing, String.split, String.compareTo, createFieldsIndexMap, BufferedReader line splitting (streamed, gzip) |

The ports take it as a path dependency (`hmf-common = { path = "../hmf-common-rs" }`), so each keeps
its own `Cargo.lock` and build folder. `cargo test --release` here runs the shared tests; each port's
regression suite is the check for a change here. Some segmentation tests read the Java tests' resource files
from `../hmf-common/src/test/resources` (or `HMF_COMMON_TEST_RESOURCES`) and are skipped without them.
Findings in this crate's own code (ids SG...) are in [REPORT.md](REPORT.md).
