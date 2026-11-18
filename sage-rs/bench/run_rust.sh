#!/usr/bin/env bash
# Run the Rust port with the same arguments as run_java.sh.
# usage: run_rust.sh <name> <threads> [extra sage args, e.g. -specific_regions "chr17:7660000-7700000"]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/env.local.sh"
name=$1; threads=$2; shift 2
# SAGE_TUMOR_ONLY=1 leaves the reference sample out
reference=(-reference "$SAGE_NORMAL" -reference_bam "$SAGE_BAM_DIR/$SAGE_NORMAL.redux.bam" -ref_sample_count 1)
if [ -n "${SAGE_TUMOR_ONLY:-}" ]; then reference=(-ref_sample_count 0); fi
bin=${SAGE_RS_BIN:-$here/../target/release/sage-rs}
out=$SAGE_OUT_DIR/rust/$name
mkdir -p "$out"
v=$SAGE_RESOURCES/dna/variants
/usr/bin/time -l "$bin" \
    -tumor "$SAGE_TUMOR" -tumor_bam "$SAGE_BAM_DIR/$SAGE_TUMOR.redux.bam" \
    "${reference[@]}" \
    -jitter_bqr_dir "$SAGE_BAM_DIR" \
    -ref_genome "$SAGE_REF_GENOME" -ref_genome_version 38 \
    -hotspots "$v/KnownHotspots.somatic.38.vcf.gz" \
    -driver_gene_panel "$SAGE_RESOURCES/common/DriverGenePanel.38.tsv" \
    -high_confidence_bed "$v/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz" \
    -ensembl_data_dir "$SAGE_RESOURCES/common/ensembl_data" \
    -include_mt \
    -threads "$threads" \
    -output_vcf "$out/$name.sage.vcf.gz" \
    "$@" > "$out/$name.log" 2>&1
tail -4 "$out/$name.log"
