#!/bin/bash
# usage: make_synthetic.sh <out dir> [seed] [unknown]: the synthetic BAMs (sorted, indexed) and pipeline files
set -eu
out=$1; seed=${2:-1}; unknown=${3:-}
python3 -I "$(dirname "$0")/make_synthetic.py" "$out" "$seed" $unknown
for s in tumor normal; do
  samtools sort -o "$out/$s.bam" "$out/$s.sam" 2>/dev/null
  samtools index "$out/$s.bam"
  rm "$out/$s.sam"
done
