#!/usr/bin/env bash
# Build the image on Google Cloud Build from the files committed at HEAD and push it.
#
#   rust-docker/cloud-build.sh <region>-docker.pkg.dev/<project>/<repository>/hmftools-rust[:tag]
#
# The tag defaults to the short commit id. Uses the active gcloud project.
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <registry path>[:tag]" >&2
    exit 2
fi

root=$(cd "$(dirname "$0")/.." && pwd)
image=$1
case "${image##*/}" in
    *:*) ;;
    *) image="$image:$(git -C "$root" rev-parse --short=10 HEAD)" ;;
esac
if [ -n "$(git -C "$root" status --porcelain -- hmf-common-rs hmf-plot-rs sage-rs redux-rs esvee-rs amber-rs cobalt-rs bam-tools-rs teal-rs lilac-rs cider-rs pave-rs purple-rs linx-rs sigs-rs virus-interpreter-rs peach-rs chord-rs neo-rs cuppa-rs qsee-rs orange-rs isofox-rs rust-docker)" ]; then
    echo "note: uncommitted changes are not part of the build; it uses HEAD" >&2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
# only what the Dockerfile copies, as committed: no build folders, nothing uncommitted
git -C "$root" archive HEAD \
    hmf-common-rs/Cargo.toml hmf-common-rs/Cargo.lock hmf-common-rs/build.rs hmf-common-rs/src hmf-common-rs/vendor hmf-common-rs/resources \
    hmf-plot-rs/Cargo.toml hmf-plot-rs/Cargo.lock hmf-plot-rs/.cargo hmf-plot-rs/src hmf-plot-rs/fonts hmf-plot-rs/resources hmf-plot-rs/examples \
    hmf-common/src/test/resources/segmentation hmf-common/src/test/resources/pcf \
    sage-rs/Cargo.toml sage-rs/Cargo.lock sage-rs/src sage-rs/resources \
    redux-rs/Cargo.toml redux-rs/Cargo.lock redux-rs/src \
    esvee-rs/Cargo.toml esvee-rs/Cargo.lock esvee-rs/build.rs esvee-rs/src esvee-rs/vendor \
    amber-rs/Cargo.toml amber-rs/Cargo.lock amber-rs/build.rs amber-rs/src amber-rs/resources amber-rs/testdata \
    cobalt-rs/Cargo.toml cobalt-rs/Cargo.lock cobalt-rs/build.rs cobalt-rs/src cobalt-rs/resources \
    cobalt/src/test/resources \
    bam-tools-rs/Cargo.toml bam-tools-rs/Cargo.lock bam-tools-rs/src bam-tools-rs/resources bam-tools-rs/testdata \
    teal-rs/Cargo.toml teal-rs/Cargo.lock teal-rs/src teal-rs/resources \
    lilac-rs/Cargo.toml lilac-rs/Cargo.lock lilac-rs/src lilac-rs/resources lilac-rs/testdata \
    cider-rs/Cargo.toml cider-rs/Cargo.lock cider-rs/src cider-rs/resources \
    pave-rs/Cargo.toml pave-rs/Cargo.lock pave-rs/build.rs pave-rs/src \
    purple-rs/Cargo.toml purple-rs/Cargo.lock purple-rs/.cargo purple-rs/build.rs purple-rs/src purple-rs/resources \
    linx-rs/Cargo.toml linx-rs/Cargo.lock linx-rs/.cargo linx-rs/src linx-rs/resources \
    sigs-rs/Cargo.toml sigs-rs/Cargo.lock sigs-rs/src \
    virus-interpreter-rs/Cargo.toml virus-interpreter-rs/Cargo.lock virus-interpreter-rs/src virus-interpreter-rs/testdata \
    peach-rs/Cargo.toml peach-rs/Cargo.lock peach-rs/src peach-rs/testdata \
    chord-rs/Cargo.toml chord-rs/Cargo.lock chord-rs/src chord-rs/resources chord-rs/testdata \
    neo-rs/Cargo.toml neo-rs/Cargo.lock neo-rs/src neo-rs/resources \
    cuppa-rs/Cargo.toml cuppa-rs/Cargo.lock cuppa-rs/.cargo cuppa-rs/src cuppa-rs/resources cuppa-rs/testdata \
    qsee-rs/Cargo.toml qsee-rs/Cargo.lock qsee-rs/.cargo qsee-rs/src qsee-rs/resources qsee-rs/testdata \
    orange-rs/Cargo.toml orange-rs/Cargo.lock orange-rs/src orange-rs/resources \
    isofox-rs/Cargo.toml isofox-rs/Cargo.lock isofox-rs/src \
    rust-docker/Dockerfile rust-docker/Dockerfile.dockerignore rust-docker/bin \
    | tar -x -C "$tmp"

echo "building $image"
gcloud builds submit "$tmp" --config "$root/rust-docker/cloudbuild.yaml" --substitutions "_IMAGE=$image"
