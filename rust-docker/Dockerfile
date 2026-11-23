# Image with the Rust ports of SAGE, REDUX and ESVEE, callable the way oncoanalyser calls the Java tools.
# Build from the repository root:
#   docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
# Commands: sage-rs, redux-rs, esvee-rs, and the wrappers sage, redux and esvee (rust-docker/bin), plus
# samtools and sambamba, which REDUX and ESVEE are given as -bamtool.

FROM rust:1-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        clang cmake libbz2-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY sage-rs sage-rs
COPY redux-rs redux-rs
COPY esvee-rs esvee-rs
RUN cd sage-rs && cargo build --release --locked \
    && cd ../redux-rs && cargo build --release --locked \
    && cd ../esvee-rs && cargo build --release --locked \
    && install -m 755 /src/sage-rs/target/release/sage-rs /src/redux-rs/target/release/redux-rs \
        /src/esvee-rs/target/release/esvee-rs /usr/local/bin/

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps samtools sambamba libbz2-1.0 liblzma5 zlib1g \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/sage-rs /usr/local/bin/redux-rs /usr/local/bin/esvee-rs /usr/local/bin/
COPY rust-docker/bin/sage rust-docker/bin/redux rust-docker/bin/esvee /usr/local/bin/
RUN sage -version && redux -version && samtools --version | head -1 && sambamba --version 2>&1 | head -2
