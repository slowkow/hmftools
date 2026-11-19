#!/usr/bin/env bash
# Convert a sample set's REDUX BAMs to CRAM so that CRAM input can be compared with BAM input.
# usage: make_cram_set.sh <out dir> [region, e.g. chr21] [samtools output options, e.g. version=3.1]
# Writes <out dir>/<sample>.redux.cram (+ .crai), links the BQR / jitter files next to them and writes
# <out dir>/env.sh, for use as SAGE_ENV=<out dir>/env.sh SAGE_BAM_EXT=cram run_{java,rust}.sh ...
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "${SAGE_ENV:-$here/env.local.sh}"
out=$1; region=${2:-}; opts=${3:-}
mkdir -p "$out"
for s in "$SAGE_TUMOR" "$SAGE_NORMAL"; do
    samtools view -@ 4 -C -T "$SAGE_REF_GENOME" ${opts:+--output-fmt-option "$opts"} \
        -o "$out/$s.redux.cram" "$SAGE_BAM_DIR/$s.redux.bam" $region
    samtools index "$out/$s.redux.cram"
    for f in "$SAGE_BAM_DIR/$s".*.tsv "$SAGE_BAM_DIR/$s".*.tsv.gz; do
        [ -e "$f" ] && ln -sf "$f" "$out/"
    done
done
{
    declare -p SAGE_TUMOR SAGE_NORMAL SAGE_REF_GENOME SAGE_RESOURCES SAGE_OUT_DIR | sed 's/^declare -x/export/;s/^declare --/export/'
    echo "export SAGE_BAM_DIR=$out"
} > "$out/env.sh"
