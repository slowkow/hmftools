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
if [ -n "$(git -C "$root" status --porcelain -- sage-rs redux-rs esvee-rs rust-docker)" ]; then
    echo "note: uncommitted changes are not part of the build; it uses HEAD" >&2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
# only what the Dockerfile copies, as committed: no build folders, nothing uncommitted
git -C "$root" archive HEAD \
    sage-rs/Cargo.toml sage-rs/Cargo.lock sage-rs/src sage-rs/resources \
    redux-rs/Cargo.toml redux-rs/Cargo.lock redux-rs/src \
    esvee-rs/Cargo.toml esvee-rs/Cargo.lock esvee-rs/build.rs esvee-rs/src esvee-rs/vendor \
    rust-docker/Dockerfile rust-docker/Dockerfile.dockerignore rust-docker/bin \
    | tar -x -C "$tmp"

echo "building $image"
gcloud builds submit "$tmp" --config "$root/rust-docker/cloudbuild.yaml" --substitutions "_IMAGE=$image"
