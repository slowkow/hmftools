# Container image for the Rust ports

One image with the Rust ports of SAGE, REDUX and ESVEE, for use in place of the Java tools' images in
oncoanalyser (https://github.com/slowkow/oncoanalyser, branch `rust-ports`, profile `rust_ports`).

Build from the repository root:

```bash
docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
```

`Dockerfile.dockerignore` limits the build context to the three crates and the wrappers.

## Contents

| Command | What it is |
|---|---|
| `sage-rs`, `redux-rs`, `esvee-rs` | the ports (`sage-rs/`, `redux-rs/`, `esvee-rs/`) |
| `sage`, `redux`, `esvee` | wrappers (`bin/`) that take the Java tools' command lines as oncoanalyser writes them |
| `samtools`, `sambamba` | given to REDUX and ESVEE as `-bamtool` |

The wrappers drop JVM options (`-Xmx...`, `-XX:...`, `-D...=...`) and turn Java main classes into
subcommands:

| Java main class | Rust command |
|---|---|
| `com.hartwig.hmftools.sage.SageApplication` (or none) | `sage-rs` |
| `com.hartwig.hmftools.sage.append.SageAppendApplication` | `sage-rs append` |
| `com.hartwig.hmftools.sage.tinc.TincAnalyser` | `sage-rs tinc` |
| `com.hartwig.hmftools.redux.ReduxApplication` (or none) | `redux-rs` |
| `com.hartwig.hmftools.esvee.EsveeApplication` (or none) | `esvee-rs` (all four steps) |
| `com.hartwig.hmftools.esvee.prep.PrepApplication` and the other three steps | `esvee-rs prep`, `assembly`, `depth`, `caller` |

Any other main class is refused with exit code 2. `-version` prints the version in Java's form, for
example `Sage version 5.1`.

`SAGE_JAVA_COMPAT=1`, `REDUX_JAVA_COMPAT=1` and `ESVEE_JAVA_COMPAT=1` (or `-java_compat`) make each port
reproduce the Java tool's output exactly; by default the ports correct the faults listed in their
REPORT.md files.
