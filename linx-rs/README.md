# linx-rs

A port of hmftools LINX 2.4 (Java) to Rust: clustering, chaining and annotation of PURPLE's structural
variants, with fusions, disruptions and the LINX driver catalog, for one sample, somatic or germline. It
takes LINX's command line and writes the same files; with `-java_compat` they equal Java's byte for
byte. Without it, the Java faults listed in [REPORT.md](REPORT.md) (LX1 to LX5) are corrected.

## Usage

```
linx-rs -sample <id> -sv_vcf <purple sv vcf> -purple_dir <dir> -ref_genome_version 38 \
    -ensembl_data_dir <dir> -known_fusion_file <csv> -driver_gene_panel <tsv> \
    -write_vis_data -write_neo_epitopes -output_dir <dir>
linx-rs -sample <id> -sv_vcf <purple germline sv vcf> -germline -ref_genome_version 38 \
    -ensembl_data_dir <dir> -driver_gene_panel <tsv> -output_dir <dir>
```

These are oncoanalyser's LINX_SOMATIC and LINX_GERMLINE steps. JVM options (`-Xmx...`, `-Xms...`,
`-XX:...`, `-D...=...`) are accepted and ignored, as oncoanalyser passes them to the bioconda wrapper.
`linx-rs -version` prints `Linx version 2.4`.

Options, as in Java:

| Option | Meaning |
|---|---|
| `-sample` | the sample id (one sample; LINX's sample lists and id files are not ported) |
| `-sv_vcf` | PURPLE's SV VCF; by default `<purple_dir>/<sample>.purple.sv.vcf.gz` (germline: `.purple.sv.germline.vcf.gz`) |
| `-purple_dir`, `-sample_data_dir` | PURPLE's output: copy number segments, purity, driver catalog and gene copy numbers (somatic) |
| `-germline` | germline mode: germline SVs with an assumed copy number of 2, germline disruptions |
| `-ref_genome_version` | 37 (default) or 38 |
| `-ensembl_data_dir` | the Ensembl data cache (required) |
| `-known_fusion_file` | known fusion data; fusions are called only with it |
| `-driver_gene_panel` | driver genes: drivers, disruptions and the driver catalog need it |
| `-write_vis_data`, `-no_vis_files`, `-write_type` | the visualiser files (written for a somatic sample unless `-no_vis_files`) |
| `-write_neo_epitopes` | the neo-epitope candidate fusions |
| `-fusion_gene_distance` | upstream distance for fusion candidate genes (default 100000) |
| `-skip_unphased_fusions`, `-write_all_vis_fusions`, `-log_reportable_fusions`, `-fusion_debug_sv_genes` | fusion options |
| `-proximity_distance` (default 5000), `-chaining_sv_limit`, `-annotations`, `-gene_id_file` | clustering and chaining options |
| `-fragile_site_file`, `-line_element_file` | in place of the bundled files (`resources/`) |
| `-output_dir`, `-log_level`, `-log_debug`, `-threads` | as in Java (LINX uses threads only for several samples) |
| `-perf_debug` | logs the time of each phase |

Port options: `-java_compat` (or `LINX_JAVA_COMPAT=1`) reproduces every Java fault;
`-java_compat_except LX2,LX5` (or `LINX_JAVA_COMPAT_EXCEPT`) all but those listed.

Not ported: the visualiser (a separate main class), RNA fusion matching (`-rna_fusions_file` is
refused), cohort and batch outputs, and the diagnostic files.

## Outputs

In `<output_dir>`, for sample `<id>` (germline: `<id>.linx.germline.*`):

| File | Content |
|---|---|
| `<id>.linx.svs.tsv` | each SV: cluster, clustering reasons, fragile site and LINE annotations, JCN range, genes, local topology |
| `<id>.linx.clusters.tsv` | each cluster: category, resolved type, SV count, description |
| `<id>.linx.links.tsv` | chain links (templated insertions) |
| `<id>.linx.breakend.tsv` | breakend transcripts (canonical, driver gene and fusion transcripts) |
| `<id>.linx.fusion.tsv` | fusions (somatic) |
| `<id>.linx.drivers.tsv`, `<id>.linx.driver.catalog.tsv` | driver events and the driver catalog (somatic) |
| `<id>.linx.germline.disruption.tsv`, `<id>.linx.germline.driver.catalog.tsv` | germline disruptions and drivers |
| `<id>.linx.neoepitope.tsv` | neo-epitope candidates (`-write_neo_epitopes`) |
| `<id>.linx.vis_*.tsv` | visualiser data: SVs, segments, copy number, gene exons, protein domains, fusions |
| `linx.version` | version and build date |

## How it works

The code follows Java's packages: `types` (SVs, breakends, links, clusters, chains), `analysis`
(clustering, classification, double minutes, cluster annotations), `chaining`, `cn` (copy number
segments, LOH and JCN), `gene` (breakend genes and transcripts), `fusion` (fusions, disruptions,
neo-epitopes), `drivers`, `vis` and `output`. Java's object graph is kept: each object sits behind an
`Rc` handle compared by identity, with its mutable fields in cells. Java's HashMaps and HashSets whose
iteration order reaches the output are reproduced: with string or integer keys by `JavaHashMap`, with
object keys (identity hashes) by `idmap.rs`, which follows a JVM with `-XX:hashCode=2` under
`-java_compat` and insertion order otherwise (finding LX1).

## Testing

```
CARGO_BUILD_JOBS=2 cargo test --release
```

runs LINX's Java unit tests, ported with their synthetic fixtures (clustering, chaining, classification,
copy number, LINE, double minutes, fusions, disruptions, drivers, germline), and tests of the corrected
faults. The comparisons with Java are described in [REPORT.md](REPORT.md).
