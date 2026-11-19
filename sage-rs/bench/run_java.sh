#!/usr/bin/env bash
# Run Java SAGE with the oncoanalyser somatic arguments (minus TINC) as the parity baseline.
# usage: run_java.sh <name> <threads> [extra sage args, e.g. -specific_regions "chr17:7660000-7700000"]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "${SAGE_ENV:-$here/env.local.sh}"   # SAGE_ENV selects another sample set
name=$1; threads=$2; shift 2
# SAGE_TUMOR_ONLY=1 leaves the reference sample out; SAGE_BAM_EXT=cram reads <sample>.redux.cram
reference=(-reference "$SAGE_NORMAL" -reference_bam "$SAGE_BAM_DIR/$SAGE_NORMAL.redux.${SAGE_BAM_EXT:-bam}" -ref_sample_count 1)
if [ -n "${SAGE_TUMOR_ONLY:-}" ]; then reference=(-ref_sample_count 0); fi
jar=${SAGE_JAR:-$here/../../sage/target/sage-5.1-jar-with-dependencies.jar}
out=$SAGE_OUT_DIR/java/$name
mkdir -p "$out"
v=$SAGE_RESOURCES/dna/variants
/usr/bin/time -l java ${JAVA_OPTS:--Xmx24g} -jar "$jar" \
    -tumor "$SAGE_TUMOR" -tumor_bam "$SAGE_BAM_DIR/$SAGE_TUMOR.redux.${SAGE_BAM_EXT:-bam}" \
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
