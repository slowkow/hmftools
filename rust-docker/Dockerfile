# Image with the Rust ports of SAGE, REDUX and ESVEE, callable the way oncoanalyser calls the Java tools.
# Build from the repository root:
#   docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
# Commands: sage-rs, redux-rs, esvee-rs, and the wrappers sage, redux and esvee (rust-docker/bin), plus
# samtools and sambamba, which REDUX and ESVEE are given as -bamtool, and procps (`ps`), which Nextflow
# uses to record each task's CPU and memory.

FROM rust:1-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        clang libclang-dev cmake libbz2-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY hmf-common-rs hmf-common-rs
# the Java test resources that the shared crate's segmentation tests read (build stage only)
COPY hmf-common/src/test/resources/segmentation hmf-common/src/test/resources/segmentation
COPY hmf-common/src/test/resources/pcf hmf-common/src/test/resources/pcf
COPY sage-rs sage-rs
COPY redux-rs redux-rs
COPY esvee-rs esvee-rs
# build each port and run its unit tests (and those of the shared crate) on this platform; a failing
# test fails the image build
RUN rustc --version \
    && (cd hmf-common-rs && cargo test --release --locked) \
    && for crate in sage-rs redux-rs esvee-rs; do \
        (cd $crate && cargo build --release --locked && cargo test --release --locked) || exit 1; \
        install -m 755 $crate/target/release/$crate /usr/local/bin/; \
    done

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps samtools sambamba libbz2-1.0 liblzma5 zlib1g \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/sage-rs /usr/local/bin/redux-rs /usr/local/bin/esvee-rs /usr/local/bin/
COPY rust-docker/bin/sage rust-docker/bin/redux rust-docker/bin/esvee /usr/local/bin/
RUN sage -version && redux -version && esvee -version && samtools --version | head -1 && sambamba --version 2>&1 | head -2
