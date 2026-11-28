# markdups-rs

Picard's MarkDuplicates in Rust, as GATK 4.6.1.0 runs it in oncoanalyser's RNA route (GATK 4.6.1.0 bundles
Picard 3.3.0 and htsjdk 4.1.3), and a mode that does `samtools sort` and MarkDuplicates in one step.

- **`markdups-rs MarkDuplicates ...`** takes GATK's command line (`--NAME value`, the short forms `-I -O -M`,
  several `--INPUT`, `--java-options` accepted) and writes the same three files: the BAM with duplicate flags and
  Picard's @PG line and `PG:Z` tags, the metrics file, and the `.bai` index (named `<output minus .bam>.bai`, as
  htsjdk names it; the pipeline renames it).
- With **`-java_compat`** (or `MARKDUPS_JAVA_COMPAT=1`) the output equals GATK's: the BAM's uncompressed bytes
  (header text, references and every record as htsjdk encodes it) and the metrics file apart from its
  "Started on" date line. The compressed bytes differ (another deflate implementation), so the index differs in
  offsets but is what htsjdk writes for this BAM (below).
- **Without it, Picard's faults are corrected** (REPORT.md, findings MD1 to MD7). `-java_compat_except MD1,MD3`
  keeps Picard's behaviour for all but the listed findings. In corrected mode the @PG line's VN says
  `markdups-rs <version> (MarkDuplicates of GATK 4.6.1.0, corrected)` instead of GATK's `Version:4.6.1.0`.
- **`markdups-rs SortAndMarkDuplicates ...`** reads STAR's unsorted BAM and writes what `samtools sort` (1.21)
  followed by MarkDuplicates writes, without the sorted intermediate BAM.

## Use

```
markdups-rs MarkDuplicates --INPUT TUMOR.bam --OUTPUT markdups_TUMOR/TUMOR.md.bam \
    --METRICS_FILE markdups_TUMOR/TUMOR.md.metrics --TMP_DIR . --CREATE_INDEX [--threads N] [-java_compat]

gatk --java-options "-Xmx30000M" MarkDuplicates --INPUT ...        # the image's wrapper, same thing

markdups-rs SortAndMarkDuplicates --INPUT Aligned.out.bam \
    --SORTED_INPUT_NAME TUMOR.bam \
    --SAMTOOLS_SORT_COMMAND_LINE "samtools sort -T TUMOR --threads 8 -o TUMOR.bam Aligned.out.bam" \
    --OUTPUT markdups_TUMOR/TUMOR.md.bam --METRICS_FILE markdups_TUMOR/TUMOR.md.metrics --TMP_DIR . --CREATE_INDEX \
    [--SORT_MEMORY_MB 2048] [--threads N] [-java_compat]
```

Threads: `--threads N`, else `MARKDUPS_THREADS`, else the CPUs the process may use (cgroup quotas included).
The output does not depend on the thread count.

### Options

All of MarkDuplicates' arguments are recognised, with GATK's defaults (COMPRESSION_LEVEL 2, as the `gatk`
launcher sets `samjdk.compression_level`; a `-Dsamjdk.compression_level=N` in `--java-options` changes it).

| Argument | Support |
|---|---|
| INPUT (several), OUTPUT, METRICS_FILE, TMP_DIR, CREATE_INDEX | yes |
| REMOVE_DUPLICATES, REMOVE_SEQUENCING_DUPLICATES, TAGGING_POLICY, CLEAR_DT | yes |
| READ_NAME_REGEX (default, `null`, or a regex with three groups: Rust `regex` syntax), OPTICAL_DUPLICATE_PIXEL_DISTANCE, MAX_OPTICAL_DUPLICATE_SET_SIZE | yes |
| DUPLICATE_SCORING_STRATEGY (SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM) | yes |
| ASSUME_SORTED, ASSUME_SORT_ORDER coordinate | yes; queryname-sorted input is refused |
| PROGRAM_RECORD_ID (or `null`), PROGRAM_GROUP_VERSION/COMMAND_LINE/NAME, ADD_PG_TAG_TO_READS, COMMENT | yes |
| VALIDATION_STRINGENCY | yes: STRICT (default) stops on the record errors htsjdk stops on, LENIENT warns, SILENT does not check |
| COMPRESSION_LEVEL | yes (deflate level of libdeflate) |
| MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, SORTING_COLLECTION_SIZE_RATIO, MAX_RECORDS_IN_RAM, VERBOSITY, QUIET, USE_JDK_DEFLATER/INFLATER, REFERENCE_SEQUENCE, the FLOW_* tuning arguments | accepted; they tune Picard's memory or logging and change no output |
| BARCODE_TAG, READ_ONE/TWO_BARCODE_TAG, MOLECULAR_IDENTIFIER_TAG, DUPLEX_UMI, TAG_DUPLICATE_SET_MEMBERS, FLOW_MODE, CREATE_MD5_FILE | refused unless at their defaults |

Not GATK's: `--threads N`; `-java_compat`, `-java_compat_except ids` (also with `--`). SortAndMarkDuplicates adds
`--SORTED_INPUT_NAME` (the name recorded as `--INPUT` in the @PG and metrics command line, the sorted BAM the
two-step pipeline would have read), `--SAMTOOLS_SORT_COMMAND_LINE` (adds samtools' @PG line with this CL, as
`samtools sort` would; without it there is none), `--SAMTOOLS_VERSION` (1.21), `--SORT_MEMORY_MB` (records held
before a sorted run is written to TMP_DIR; default 2048).

Input: BAM (BGZF; STAR's uncompressed `--outBAMcompression 0` BAM included). SAM and CRAM are refused.

## How it works

Two passes, like Picard. Pass 1 reads the input and builds a fragment end (library, 5' unclipped position,
strand, score, tile/x/y from the read name) per primary mapped read and a pair end once both mates are seen;
Picard sorts all of them and walks the sorted lists in chunks of equal position. Every chunk holds ends of one
reference, so the port sorts and walks each reference's ends when the input moves past it (pairs with mates on two
references at the end), which keeps memory to one reference's ends. Duplicates are a bit per record. Pass 2
reads the input again and writes it with the flags, the DT and PG tags, counting the metrics, deflating 64 KiB
blocks on worker threads and building the index from the blocks' addresses.

Records are written as htsjdk writes them when an attribute was set (Picard clears DT and sets PG on every
record): optional fields sorted by htsjdk's tag key (second character first) and integers in the smallest type,
the bin recomputed, missing qualities as 0xff, the unused half of an odd last base byte zeroed. The header is
re-encoded as htsjdk does (@HD VN:1.6 first; @SQ, @RG, @PG, @CO grouped).

The index follows htsjdk's BAMIndexer as BAMFileWriter drives it with CREATE_INDEX: the same bins, chunk merging,
metadata pseudo-bin and linear index for this BAM. It equals htsjdk's BuildBamIndex of the port's BAM except for the
last record's chunk end, which the writer gives inside the last data block and BuildBamIndex as the start of the
end-of-file block (GATK's own index of its own BAM shows the same). htslib (samtools) reads it.

SortAndMarkDuplicates sorts as samtools 1.21 does (reference index unsigned, so unplaced reads last; position;
reverse strand; input order among equals), changes the header as samtools does (SO:coordinate, GO and SS removed,
its @PG line per chain end), keeps the last run in memory and writes earlier runs to TMP_DIR at deflate level 1,
then runs both passes over the merged runs. It saves the sorted BAM (written at samtools' default level and read
twice by GATK) and its disk space.

## Outputs, identical under -java_compat

Checked by the comparison scripts of the port's work folder (REPORT.md has the results): records, header and uncompressed bytes against
GATK 4.6.1.0 (the pipeline's image, and the same jar on Java 23), at 1 and 2 threads, on a tumour RNA slice, three
public HG008 slices, 39 synthetic cases (one per rule, plus several inputs), and the one-step mode against
samtools 1.21 sort + GATK.

## Tests

`cargo test --release`: unit tests of each module and Picard's tests of the ported classes
(`src/picard_tests.rs`: AbstractMarkDuplicatesCommandLineProgramTest's and MarkDuplicatesTest's cases through
a tester like Picard's SamFileTester; OpticalDuplicateFinderTest and ReadNameParserTest cases in
`src/optical.rs`; DuplicationMetrics' library sizes against Picard's own `DuplicationMetrics.main`).

## What comes from Picard and htsjdk

Ported from Picard 3.3.0 (MIT licence, `LICENSE-picard`; Copyright (c) 2017-2021 Broad Institute):
MarkDuplicates (`markdups.rs`), AbstractMarkDuplicatesCommandLineProgram's metrics and optical duplicate tracking,
ReadEnds / ReadEndsForMarkDuplicates and their comparator, LibraryIdGenerator, OpticalDuplicateFinder,
ReadNameParser and GraphUtils (`optical.rs`), DuplicationMetrics (`metrics.rs`), the argument definitions
(`args.rs`), and the tests in `picard_tests.rs`. Ported from htsjdk 4.1.3 (MIT licence, `LICENSE-htsjdk`): BAMRecordCodec,
BinaryTagCodec, SAMTextHeaderCodec, SAMRecord.isValid and the CIGAR checks (`validate.rs`), BAMIndexer,
BinningIndexBuilder and BinaryBAMIndexWriter (`bai.rs`), MetricsFile and FormatUtil, SamFileHeaderMerger and
MergingSamRecordIterator (`input.rs`), DuplicateScoringStrategy and Murmur3. The samtools sort order and header
changes (`sort.rs`) follow samtools 1.21 (MIT/Expat). The rest of the crate is GPL-3.0, as the repository.

## Building

`CARGO_BUILD_JOBS=2 cargo build --release`; depends on `../hmf-common-rs` (path dependency) and links libdeflate
through `libdeflate-sys` (the copy rust-htslib uses). The image needs `markdups-rs` on the PATH and
`rust-docker/bin/gatk` as `gatk`.

## Upstream

Ports Picard 3.3.0 with htsjdk 4.1.3, as bundled in GATK 4.6.1.0 (the image oncoanalyser runs); `-java_compat`
reproduces that GATK. Later Picard releases are followed as hmftools is, change by change (see
[`../UPSTREAM.md`](../UPSTREAM.md)).
