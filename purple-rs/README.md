# purple-rs

PURPLE in Rust: a port of hmftools PURPLE 4.5 (`com.hartwig.hmftools.purple.PurpleApplication`), the
purity, ploidy and copy-number step of oncoanalyser. It reads the AMBER, COBALT, ESVEE and PAVE outputs
and writes the same files as the Java tool. With `-java_compat` the output is byte-identical to Java's
(the VCFs, their BGZF bytes and tabix indexes too); by default the faults found in the Java code are corrected
([REPORT.md](REPORT.md), PU1-PU9).

## Running

```
purple-rs <the Java options>
```

JVM options (`-Xmx...`, `-Xms...`, `-XX:...`) are ignored. `purple-rs -version` prints
`Purple version 4.5`; `-help` lists the options.

oncoanalyser's tumor/normal call:

```
purple-rs -tumor T -reference N -amber amber/ -cobalt cobalt/ -esvee_dir esvee/ \
  -pave_somatic_dir pave/somatic/ -pave_germline_dir pave/germline/ \
  -ref_genome ref.fasta -ref_genome_version 38 -driver_gene_panel DriverGenePanel.38.tsv \
  -ensembl_data_dir ensembl_data/ -somatic_hotspots KnownHotspots.somatic.38.vcf.gz \
  -germline_hotspots KnownHotspots.germline.38.vcf.gz \
  -germline_amp_del_freq_file cohort_germline_amp_del_freq.38.csv \
  -circos circos -threads 2 -output_dir purple/
```

Modes, as in Java: tumor/normal (`-tumor` and `-reference`), tumor-only (`-tumor` alone), germline-only
(`-reference` alone; the inputs are named after the reference), targeted (any of these with
`-target_regions_bed`).

## Options

| Option | Default | Meaning |
|---|---|---|
| `-tumor`, `-reference` | | sample ids (at least one) |
| `-output_dir` | | output directory (under `-sample_dir` when that is given) |
| `-amber_dir` / `-amber`, `-cobalt_dir` / `-cobalt` | | AMBER and COBALT outputs |
| `-esvee_dir`, `-pave_somatic_dir`, `-pave_germline_dir` | | ESVEE and PAVE outputs (`<sample>.esvee.somatic.vcf.gz`, ...) |
| `-somatic_sv_vcf`, `-germline_sv_vcf`, `-somatic_vcf`, `-germline_vcf` | | the VCFs named directly |
| `-sample_dir` | | pipeline directory holding the tools' outputs (OA_V3_1 layout) |
| `-redux_tumor_dir` | | REDUX outputs (the MSI prediction used in targeted mode) |
| `-ref_genome`, `-ref_genome_version` | V37 | reference FASTA and version (37 or 38) |
| `-driver_gene_panel`, `-ensembl_data_dir` | | driver genes; Ensembl cache (required) |
| `-somatic_hotspots`, `-germline_hotspots` | | hotspot VCFs |
| `-germline_amp_del_freq_file` | | cohort frequencies of germline amplifications and deletions |
| `-target_regions_bed`, `-target_coding_factor`, `-target_tmb_ratio`, `-target_tml_ratio` | | targeted mode |
| `-min_purity`, `-max_purity`, `-purity_increment`, `-min_ploidy`, `-max_ploidy`, `-ploidy_penalty_*`, `-deviation_penalty_*`, `-ambiguous_baf_limit`, ... | Java's | the fitting parameters (see `-help`) |
| `-somatic_min_peak`, `-somatic_min_variants`, `-somatic_min_purity`, `-somatic_penalty_weight`, `-highly_diploid_percentage` | Java's | somatic fit parameters |
| `-sv_qual_filter` | 0 | minimum SV QUAL (0: off) |
| `-write_all_somatics` | | in tumor-only mode, write filtered somatic variants too |
| `-circos` | | as in Java, turns the Circos plots on: the Circos input files are written and `circos.png` and `input.png` drawn natively (the binary named is not run) |
| `-circos_native` | | the same without naming a binary (`-circos` must name an existing file, as in Java) |
| `-external_circos` | | run the Circos binary given by `-circos` instead, as Java does |
| `-no_charts`, `-ignore_plot_errors`, `-skip_circos_genes` | | charting switches, as in Java |
| `-write_plot_data` | | keep `plot/<tumor>.somatic_data.tsv`, the rainfall plot's data (Java writes it for R and deletes it) |
| `-specific_chr`, `-specific_regions` | | restrict to chromosomes or regions |
| `-threads` | 1 | worker threads |
| `-log_level`, `-log_debug` | INFO | log level |
| `-java_compat`, `-java_compat_except PU1,...` | | Java's behaviour for every finding, or all but those listed; also `PURPLE_JAVA_COMPAT=1`, `PURPLE_JAVA_COMPAT_EXCEPT` |

## Outputs

In `-output_dir`, `<tumor>.purple.<name>`:

| File | Content |
|---|---|
| `purity.tsv`, `qc` | the fit (purity, ploidy, method, MSI/TMB/TML), QC status |
| `purity.range.tsv` | every candidate fit by score |
| `segment.tsv`, `pre_resegmented_regions.tsv` | fitted regions; regions before resegmentation |
| `cnv.somatic.tsv`, `cnv.gene.tsv`, `chromosome_arm.tsv` | copy-number segments, genes, chromosome arms |
| `somatic.vcf.gz`, `germline.vcf.gz`, `sv.vcf.gz`, `sv.germline.vcf.gz` | the annotated VCFs (tabix indexed) |
| `driver.catalog.somatic.tsv`, `driver.catalog.germline.tsv` | driver catalogs |
| `germline_amp_del.tsv` | germline amplifications and deletions |
| `somatic.clonality.tsv`, `somatic.hist.tsv` | peak model; variant copy number histogram |
| `target_region_cn.tsv` | targeted mode: copy number of each COBALT window over a target region |
| `purple.version` | version and build date |
| `circos/` | with `-circos`: the Circos data and configuration files, gaps.txt |
| `plot/` | unless `-no_charts` or germline-only: `copynumber.png` and `map.png` (when there are autosomal segments with BAF support), `purity.range.png`, `segment.png`; with a somatic VCF also `somatic.png`, `somatic.clonality.png`, `somatic.rainfall.png`; with `-circos` also `circos.png` and `input.png` |

The plots are drawn natively by [hmf-plot-rs](../hmf-plot-rs/README.md) (BUI's software raster), without R,
Perl or Circos, from the same files the R scripts and Circos read, with the same sizes, layout, scales,
colours and labels. `src/plots.rs` holds the plot-specific steps (the R scripts' dplyr and ggplot2 stats).

## Corrected by default

See [REPORT.md](REPORT.md) for each finding. In short: kataegis ids numbered in genome order whatever the
thread count (PU1); the AMBER points of baf.circos in genome order (PU2); somatic VCFs checked for the
tumor sample in tumor-only mode (PU4); arm maxima of all-negative arms (PU6); the phased-inframe
reportability check (PU8); germline amplification drivers typed GERMLINE_AMP (PU9). In the plots: every SNV
point drawn on the Circos plot (PL1), the bar outlines the R scripts ask for (PL2), rainfall distances within
a chromosome (PL3), and up to 100,000 rainfall SNVs (PL4).

## Building and testing

```
CARGO_BUILD_JOBS=2 cargo build --release
CARGO_BUILD_JOBS=2 cargo test --release
```

Depends on `../hmf-common-rs` and `../hmf-plot-rs` (path dependencies; hmf-plot-rs builds BUI's C kit, which
needs `cmake`). The unit tests are ports of the Java tests with
synthetic data.

## Upstream

Started from hmftools commit `ca04e6c39d` (2026-11-20); `-java_compat` reproduces the Java built at that commit. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
