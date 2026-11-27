# amber-rs

A Rust port of hmftools AMBER 4.4 (`../amber`): allele counts at known heterozygous germline sites in a
tumor and/or reference BAM or CRAM, the B-allele frequencies (BAF) and their segmentation that PURPLE
reads, regions of homozygosity, the tumor's contamination and, in tumor-only mode, the purity analysis
and LOH recovery. It takes the same options and writes the same files as the Java `AmberApplication`.
With `-java_compat` the files are byte-identical to Java's, the gzip/BGZF bytes and the tabix index too.
Faults found in the Java code are corrected by default (findings AM1-AM7 and AM51-AM55 in
[REPORT.md](REPORT.md)).

The oncoanalyser fork runs Hartwig's `amber:4.4-beta.1` image, built from the tag `amber-v4.4-beta.1`,
whose AMBER sources are the ones this port follows (upstream's module files name AMBER 4.3, which lacks
the BAF file's `lohRecovered` column; PURPLE reads that file by column name).

## Build and run

```
cargo build --release          # target/release/amber-rs
amber-rs -tumor TUMOR -tumor_bam tumor.bam -reference NORMAL -reference_bam normal.bam \
    -loci AmberGermlineSites.38.tsv.gz -ref_genome_version 38 -threads 8 -output_dir amber/
```

Modes:
- tumor/normal;
- tumor-only (no reference: the purity analysis and LOH recovery);
- germline-only (no tumor);
- targeted (`-target_regions_bed`).

Several reference BAMs can be given (`-reference A,B -reference_bam a.bam,b.bam`). JVM options are
accepted and ignored, so the bioconda wrapper's command line works. `amber-rs -version` prints
`Amber version 4.4`.

## Options

As Java's (`amber-rs -help` lists them with their defaults):
- **Inputs:** `-tumor`, `-tumor_bam`, `-reference`, `-reference_bam`, `-loci` (TSV or VCF),
  `-ref_genome` (required for CRAM), `-ref_genome_version`.
- **Filters:**
  - `-min_base_quality`, `-min_map_quality`;
  - `-min_depth_percent`, `-max_depth_percent`;
  - `-min_het_af_percent`, `-max_het_af_percent`;
  - `-tumor_min_depth`, `-tumor_only_min_support`, `-tumor_only_min_vaf`.
- **Restrictions and reading:** `-specific_chr`, `-specific_regions`, `-target_regions_bed`,
  `-position_gap`, `-sequencing_type` (ILLUMINA, ULTIMA, SBX).
- **Output:** `-write_tumor_data`, `-write_unfiltered_germline`, `-skip_baf_segmentation`, `-output_dir`.
- **Running:** `-threads` (0 is taken as 1; Java fails, AM7), `-log_debug`, `-log_level`,
  `-bam_validation` (accepted, not applied).
- **Java behaviour:** `-java_compat` gives Java's behaviour for every finding;
  `-java_compat_except AM5,...` gives it for all but those listed. The environment variables
  `AMBER_JAVA_COMPAT=1` and `AMBER_JAVA_COMPAT_EXCEPT` do the same.

## Outputs

In `-output_dir`, as Java writes them:
- `<tumor>.amber.baf.tsv.gz`, with the BAF segments in `<tumor>.amber.baf.pcf`;
- `<tumor>.amber.qc` and `<tumor>.amber.contamination.tsv`;
- `<reference>.amber.homozygousregion.tsv` and `<reference>.amber.snp.vcf.gz` with its `.tbi`;
- `<tumor>.amber.tumor.raw.tsv.gz` (with `-write_tumor_data`);
- `amber.version`, with the Rust build's time as `build.date`; set `AMBER_BUILD_DATE` at build time to
  change it.

## Tests

`cargo test --release` runs the ported Java unit tests, the Java end-to-end scenarios on the Java tests'
synthetic BAMs (`testdata/e2e`), and the port's own tests. Five ignored tests compare against files written
by Java and run when the files are given:
- `poisson::tests::poisson_matches_java` (`AMBER_POISSON_REF`);
- `contamination::tests::contamination_matches_java` (`AMBER_CONTAMINATION_REF`);
- `purity::binomial::tests::java_grid` (`AMBER_BINOMIAL_GRID`);
- `purity::hash_order::tests::java_orders` (`AMBER_HASH_ORDERS`);
- `purity::driver::tumor_only_driver`, the Rust side of the step-by-step tumor-only comparison.

The comparisons with Java are in [REPORT.md](REPORT.md).
