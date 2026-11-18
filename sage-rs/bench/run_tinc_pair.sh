#!/usr/bin/env bash
# Run the standalone TINC analyser in Java and Rust on copies of the same SAGE VCF and compare the rewritten VCFs.
# usage: run_tinc_pair.sh <name> <input.sage.vcf.gz> [threads]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "${SAGE_ENV:-$here/env.local.sh}"   # SAGE_ENV selects another sample set
name=$1; input=$2; threads=${3:-4}
jar=${SAGE_JAR:-$here/../../sage/target/sage-5.1-jar-with-dependencies.jar}
v=$SAGE_RESOURCES/dna/variants
for t in java rust; do
    d=$SAGE_OUT_DIR/$t/$name
    mkdir -p "$d"
    cp -f "$input" "$d/in.sage.vcf.gz"
    cp -f "$input.tbi" "$d/in.sage.vcf.gz.tbi"
    args=(-input_vcf "$d/in.sage.vcf.gz" -tumor "$SAGE_TUMOR" -reference "$SAGE_NORMAL" -ref_genome_version 38
        -pon_file "$v/hmf_wgs_sage_pon_1000.38.tsv.gz" -gnomad_freq_dir "$v/gnomad/"
        -write_tinc_vcf -write_fit_variants -threads "$threads")
    if [ $t = java ]; then
        /usr/bin/time -l java -Xmx16g -cp "$jar" com.hartwig.hmftools.sage.tinc.TincAnalyser "${args[@]}" \
            -ref_genome "$SAGE_REF_GENOME" -log_debug > "$d/log" 2>&1
    else
        /usr/bin/time -l "$here/../target/release/sage-rs" tinc "${args[@]}" > "$d/log" 2>&1
    fi
    echo "== $t"
    grep -i "tinc level\|fitCount\|recovered\| real \|peak memory\|rror" "$d/log" | sed 's/^[0-9:.]* \[[A-Z ]*\] //' | cut -c1-200
done
j=$SAGE_OUT_DIR/java/$name; r=$SAGE_OUT_DIR/rust/$name
cmp <(gzcat "$j/in.sage.tinc.vcf.gz") <(gzcat "$r/in.sage.tinc.vcf.gz") && echo "tinc VCF: IDENTICAL"
cmp "$j/in.sage.fit_variant.tsv" "$r/in.sage.fit_variant.tsv" && echo "fit variants: IDENTICAL"
