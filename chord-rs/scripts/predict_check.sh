#!/bin/bash
# The prediction alone, R against the port: predict_check.sh <rows> <seed> <dir>. Runs chord_predict.R
# (as Java runs it) and `chord-rs predict -java_compat` with the embedded forest and with CHORD.rds read
# directly, on a table of random synthetic samples, and compares the prediction files. R needs the
# randomForest package; for rows whose BRCA1 and BRCA2 probabilities tie, R picks at random (CH1), and the
# defined reference is R after set.seed(1): give R_PROFILE_USER a profile with `set.seed(1)` for that.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
ROWS=$1; SEED=$2; D=$3
mkdir -p "$D"
python3 -I "$here/random_contexts.py" "$ROWS" "$SEED" "$D/contexts.tsv"
RDS=$root/chord/src/main/R/CHORD.rds
BIN=${CHORD_RS_BIN:-$root/chord-rs/target/release/chord-rs}
Rscript "$root/chord/src/main/R/chord_predict.R" "$RDS" "$D/contexts.tsv" "$D/r.tsv" INFO > "$D/r.log" 2>&1
"$BIN" predict -mut_contexts_file "$D/contexts.tsv" -output_file "$D/rust.tsv" -java_compat > "$D/rust.log" 2>&1
"$BIN" predict -mut_contexts_file "$D/contexts.tsv" -output_file "$D/rust_rds.tsv" -chord_model_file "$RDS" -java_compat > "$D/rust_rds.log" 2>&1
"$BIN" predict -mut_contexts_file "$D/contexts.tsv" -output_file "$D/rust_fixed.tsv" > "$D/rust_fixed.log" 2>&1
cmp "$D/r.tsv" "$D/rust.tsv" && echo "embedded forest: identical ($ROWS rows)"
cmp "$D/r.tsv" "$D/rust_rds.tsv" && echo "CHORD.rds read directly: identical"
echo "corrected mode: $(diff "$D/r.tsv" "$D/rust_fixed.tsv" | grep -c '^>') of $ROWS rows differ"
cut -f5,6 "$D/r.tsv" | sort | uniq -c
