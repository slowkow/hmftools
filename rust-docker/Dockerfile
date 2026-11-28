# Image with the Rust ports of SAGE, REDUX, ESVEE, AMBER, COBALT, BamTools (metrics), TEAL, LILAC,
# CIDER, PAVE, PURPLE, LINX, SIGS, Virus Interpreter, PEACH, CHORD, NEO, CUPPA, qsee, ORANGE and ISOFOX, and of Picard's MarkDuplicates (as GATK 4.6.1.0 runs it), callable the way oncoanalyser calls the Java tools.
# Build from the repository root:
#   docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
# Commands: sage-rs, redux-rs, esvee-rs, amber-rs, cobalt-rs, bam-tools-rs, teal-rs, lilac-rs, cider-rs, pave-rs,
# purple-rs, linx-rs, sigs-rs, virus-interpreter-rs, peach-rs, chord-rs, neo-rs, cuppa-rs, qsee-rs, orange-rs, isofox-rs, markdups-rs, and the wrappers sage, redux, esvee, amber,
# cobalt, bamtools, teal, lilac, cider, pave, purple, linx, sigs, virusinterpreter, peach, chord, neo, cuppa, qsee, orange, isofox and gatk (rust-docker/bin; gatk runs only MarkDuplicates), a `python`
# that only runs `python -m cuppa.predict` (by cuppa-rs), plus
# samtools and sambamba, which REDUX and ESVEE are given as -bamtool, and procps (`ps`), which Nextflow
# uses to record each task's CPU and memory.

# The toolchain is pinned to the one every port was checked with (rustc 1.97.1). rust:1-bookworm moved to rustc
# 1.100.0, with which a release build of cuppa-rs typed every SNV of its test VCF as UNDEFINED (the code is safe;
# a debug print in the function made it right again): a compiler fault or undefined behaviour elsewhere, not yet
# found (cuppa-rs REPORT.md, open points).
FROM rust:1.97-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        clang libclang-dev cmake libbz2-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY hmf-common-rs hmf-common-rs
# the drawing crate (BUI's raster, fetched from GitHub at a pinned commit; cmake builds its C kit) for purple-rs's plots
COPY hmf-plot-rs hmf-plot-rs
# the Java test resources that the shared crate's segmentation tests read (build stage only)
COPY hmf-common/src/test/resources/segmentation hmf-common/src/test/resources/segmentation
COPY hmf-common/src/test/resources/pcf hmf-common/src/test/resources/pcf
# the synthetic BAMs of the Java COBALT tests, which cobalt-rs's end-to-end tests read (build stage only)
COPY cobalt/src/test/resources cobalt/src/test/resources
COPY sage-rs sage-rs
COPY redux-rs redux-rs
COPY esvee-rs esvee-rs
COPY amber-rs amber-rs
COPY cobalt-rs cobalt-rs
COPY bam-tools-rs bam-tools-rs
COPY teal-rs teal-rs
COPY lilac-rs lilac-rs
COPY cider-rs cider-rs
COPY pave-rs pave-rs
COPY purple-rs purple-rs
COPY linx-rs linx-rs
COPY sigs-rs sigs-rs
COPY virus-interpreter-rs virus-interpreter-rs
COPY peach-rs peach-rs
COPY chord-rs chord-rs
COPY neo-rs neo-rs
COPY cuppa-rs cuppa-rs
COPY qsee-rs qsee-rs
COPY orange-rs orange-rs
COPY isofox-rs isofox-rs
COPY markdups-rs markdups-rs
# build each port and run its unit tests (and those of the shared crate) on this platform; a failing
# test fails the image build
RUN rustc --version \
    && (cd hmf-common-rs && cargo test --release --locked && cargo test --release --locked --features bwa) \
    && (cd hmf-plot-rs && cargo test --release --locked --all-features) \
    && for crate in sage-rs redux-rs esvee-rs amber-rs cobalt-rs bam-tools-rs teal-rs lilac-rs cider-rs pave-rs purple-rs linx-rs \
            sigs-rs virus-interpreter-rs peach-rs chord-rs neo-rs cuppa-rs qsee-rs orange-rs isofox-rs markdups-rs; do \
        (cd $crate && cargo build --release --locked && cargo test --release --locked) || exit 1; \
        install -m 755 $crate/target/release/$crate /usr/local/bin/; \
    done

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps samtools sambamba libbz2-1.0 liblzma5 zlib1g \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /usr/local/bin/sage-rs /usr/local/bin/redux-rs /usr/local/bin/esvee-rs /usr/local/bin/amber-rs \
    /usr/local/bin/cobalt-rs /usr/local/bin/bam-tools-rs /usr/local/bin/teal-rs \
    /usr/local/bin/lilac-rs /usr/local/bin/cider-rs /usr/local/bin/pave-rs \
    /usr/local/bin/purple-rs /usr/local/bin/linx-rs /usr/local/bin/sigs-rs /usr/local/bin/virus-interpreter-rs \
    /usr/local/bin/peach-rs /usr/local/bin/chord-rs /usr/local/bin/neo-rs \
    /usr/local/bin/cuppa-rs /usr/local/bin/qsee-rs /usr/local/bin/orange-rs /usr/local/bin/isofox-rs /usr/local/bin/markdups-rs /usr/local/bin/
COPY rust-docker/bin/sage rust-docker/bin/redux rust-docker/bin/esvee rust-docker/bin/amber rust-docker/bin/cobalt \
    rust-docker/bin/bamtools rust-docker/bin/teal rust-docker/bin/lilac \
    rust-docker/bin/cider rust-docker/bin/pave rust-docker/bin/purple rust-docker/bin/linx \
    rust-docker/bin/sigs rust-docker/bin/virusinterpreter rust-docker/bin/peach rust-docker/bin/chord rust-docker/bin/neo \
    rust-docker/bin/cuppa rust-docker/bin/python rust-docker/bin/qsee rust-docker/bin/orange rust-docker/bin/isofox rust-docker/bin/gatk /usr/local/bin/
RUN sage -version && redux -version && esvee -version && amber -version && cobalt -version && bamtools -version && teal -version && lilac -version && cider -version && pave -version && purple -version && linx -version && sigs -version && virusinterpreter -version && peach -version && chord -version && neo -version && cuppa -version && qsee -version && orange -version && isofox -version && gatk --version && samtools --version | head -1 && sambamba --version 2>&1 | head -2
