# hmf-common-rs

Code shared by the Rust ports of hmftools (`sage-rs`, `redux-rs`, `esvee-rs`, `amber-rs`, `cobalt-rs`,
`bam-tools-rs`, `teal-rs`, `lilac-rs`, `cider-rs`, `pave-rs`, `purple-rs`, `linx-rs`, `sigs-rs`, `virus-interpreter-rs`, `peach-rs`, `chord-rs`, and the ports that follow).
It plays the part of the Java `hmf-common` module, plus the Java library behaviour the ports reproduce
so that their output equals the Java tools'.

| Module | What it is |
|---|---|
| `args` | command lines in the style of the Java ConfigBuilder (`-name value`, `-flag`) |
| `compat` | the switch between a port's corrected behaviour and Java's (`-java_compat`, `-java_compat_except`, `<TOOL>_JAVA_COMPAT`); each port keeps its own finding ids |
| `faidx`, `refgenome` | indexed FASTA access, 1-based inclusive positions as the Java RefGenomeInterface |
| `fmt` | numbers as Java's `String.format`, `Double.toString` and `DecimalFormat` write them (`decimal_format` ports the JDK's FloatingDecimal digits and DigitList rounding; `decimal_format_fixed` for patterns such as "0.0000") |
| `log` | log lines in the layout of the Java tools; the `info!`, `warn!`, `error!` and `debug!` macros |
| `regions` | chromosome ranks and names, regions (with BaseRegion's positionWithin and positionsWithin), `-specific_regions` / `-specific_chr`; RefGenomeVersion, BaseRegion, region and BED file loading with Java's merging, PartitionUtils |
| `java_hash`, `java_id_hash` | java.util.HashMap and HashSet iteration order (with `retain` for retainAll / removeIf, which keep the table size, and `compute_if_absent`, which puts a new key at the head of its bucket as Collectors.groupingBy does), for hash codes and for identity hashes under `-XX:hashCode=2` |
| `java_sort`, `java_math` | TimSort as `List.sort`, and Java numeric behaviour such as `DoubleStream.sum`, `Math.max`, `Math.round`, `Double.compare` |
| `java_fastmath` | commons-math3 `FastMath.log` and `FastMath.exp`, which differ from `Math.log` and `Math.exp` in the last bit for some inputs (COBALT takes logs with it; AMBER's Poisson and binomial distributions take both) |
| `segmentation` | hmf-common `segmentation`, `segmentation.copynumber` and `utils.pcf`: the per-arm piecewise constant fit (a port of the R copynumber package's pcf), its penalty, the segments file (AMBER absolute, COBALT windowed) and the .pcf readers (PURPLE) |
| `ensembl` | the Ensembl data cache: genes, transcripts and exons as EnsemblDataCache loads them (canonical only, restricted gene ids, synonyms), the name and id maps in Java HashMap order; protein domains and splice acceptors, findGeneRegions, findGenesByRegion, getProteinDomainPositions, getExonRankings, the biotype constants |
| `driver_panel` | the driver gene panel file (read, write) and DriverGeneRegions (the coding exon regions of the panel's canonical transcripts) |
| `metrics` | the BamTools metrics files: summary, coverage, fragment lengths, flag counts (BamFlagStats), gene and exon depth (GeneDepthFile); writers and readers |
| `java_text` | Java's Integer/Long/Double/Boolean parsing, String.split, String.compareTo, createFieldsIndexMap, BufferedReader line splitting (streamed, gzip) |
| `bwa` (feature `bwa`) | BWA-MEM alignment as gatk-bwamem-jni does it (ESVEE, CIDER): an index image mapped from a file (read through first on a network file system, or always if asked) or held in memory, the `mem_opt_t` of `mem_opt_init()` as a Rust struct to change field by field (as the Java aligners do through a ByteBuffer), and the BwaMemAlignment records. The bwa sources (0.7.16a as in gatk-bwamem-jni) are in `vendor/bwa` and built by `build.rs` only with the feature |
| `codon` | Codons (amino acids from DNA, '.' for an unknown codon) and Nucleotides / htsjdk SequenceUtil complements and reverse complements |
| `blosum` | BlosumMapping: BLOSUM62 with the stop codon penalty |
| `csv` | commons-csv's CSVPrinter as the Kotlin tools write TSV: MINIMAL quoting, null string |
| `kotlin_random` | kotlin.random.Random(seed), the XorWow generator (nextInt, nextBits, nextDouble) |
| `sam_header` | SAM header lines as htsjdk writes them: `@HD VN:1.6` (SAMFileHeader.CURRENT_VERSION) whatever the version of the header read, as REDUX, TEAL and CIDER write their BAMs |
| `sv` | structural variants from a VCF: StructuralVariantType, the record reader with htsjdk's attribute rules, StructuralVariantFactory, EnrichedStructuralVariantFactory, StructuralVariantData (LINX) |
| `purple_files` | PURPLE's files: PurpleCopyNumber, PurityContext (purity, ploidy, gender, WGD), GeneCopyNumber, DriverCatalog (read and write), `write_lines` as Files.write |
| `known_fusion` | KnownFusionCache and KnownFusionData (the known fusion file), KnownFusionType, FusionCommon's constants |
| `transcript_utils` | TranscriptUtils: coding bases and phases |
| `poisson` | re-exports `commons_math`'s PoissonDistribution, logGamma and regularized gamma functions (AMBER, LINX) |
| `gz` | fast gzip / BGZF text reading (zlib-rs, one member after another, optionally inflating on its own thread) and a line reader (PAVE's resource files) |
| `vcf` | VCF headers, records and the VCFEncoder as htsjdk 2.23 has them, and `VcfFileWriter`: the BGZF stream (zlib level 5, as java.util.zip) with htsjdk's on-the-fly tabix index (PAVE) |
| `variant_context` | VCF records with htsjdk's typed attributes and lazy genotypes (VariantContext, LazyGenotypesContext, VCFEncoder.formatVCFField), written through `vcf`'s BGZF stream and tabix index (PURPLE) |
| `commons_math` | commons-math3 3.6.1: Gamma, Beta, Binomial, Poisson and Normal distributions, Erf.erfc, StandardDeviation, bit for bit (PURPLE; `poisson` re-exports its Poisson for AMBER and LINX) |
| `java_strictmath` | StrictMath exp, log, log10 and pow (FdLibm) |
| `kernel_estimator` | hmf-common KernelEstimator |
| `config_builder` | hmf-common ConfigBuilder: registered items, hasValue, path and number checks, -version and -help (PURPLE) |
| `ref_coords` | RefGenomeCoordinates: the GRCh37 and GRCh38 lengths of chromosomes 1-22, X and Y (resources refgenome/lengths.*.tsv) |

The ports take it as a path dependency (`hmf-common = { path = "../hmf-common-rs" }`, with
`features = ["bwa"]` for the aligner), so each keeps
its own `Cargo.lock` and build folder. `cargo test --release` here runs the shared tests; each port's
regression suite is the check for a change here. Some segmentation tests read the Java tests' resource files
from `../hmf-common/src/test/resources` (or `HMF_COMMON_TEST_RESOURCES`) and are skipped without them.
Findings in this crate's own code (ids SG...) are in [REPORT.md](REPORT.md).
