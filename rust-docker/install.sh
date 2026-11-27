#!/usr/bin/env bash
# Build the Rust ports and put them, with the sage/redux/esvee/amber/cobalt/bamtools/teal wrapper commands, into
# one folder, for oncoanalyser's --rust_ports_bin (runs these steps on the host instead of in a container).
#
#   rust-docker/install.sh <folder>
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <folder>" >&2
    exit 2
fi

root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$1"
dest=$(cd "$1" && pwd)

for crate in sage-rs redux-rs esvee-rs amber-rs cobalt-rs bam-tools-rs teal-rs; do
    (cd "$root/$crate" && cargo build --release --locked)
    install -m 755 "$root/$crate/target/release/$crate" "$dest/"
done
for tool in sage redux esvee amber cobalt bamtools teal; do
    install -m 755 "$root/rust-docker/bin/$tool" "$dest/"
done

PATH="$dest:$PATH"
sage -version
redux -version
amber -version
cobalt -version
bamtools -version
teal -version
for tool in samtools sambamba; do
    # REDUX and ESVEE are given these as -bamtool; the ports only check that the path exists
    command -v "$tool" > /dev/null || echo "warning: $tool is not on PATH; oncoanalyser passes it to REDUX and ESVEE as -bamtool" >&2
done
echo "installed in $dest; run oncoanalyser with --rust_ports_bin $dest"
