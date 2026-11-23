# Running the Rust ports in oncoanalyser

The Rust ports of SAGE, REDUX and ESVEE, packaged for use in place of the Java tools in oncoanalyser
(https://github.com/slowkow/oncoanalyser, branch `rust-ports`, profile `rust_ports`). Two forms:

**A folder on the host**, for running these steps without a container (oncoanalyser's
`--rust_ports_bin`). The ports need nothing beyond the operating system. This builds them and copies
them with the wrappers into one folder:

```bash
rust-docker/install.sh ~/rust-ports
```

**A container image**, for executors that run every task in a container (Google Cloud Batch, AWS
Batch, Kubernetes) or when Docker is used anyway. Build from the repository root:

```bash
docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
```

`Dockerfile.dockerignore` limits the build context to the three crates and the wrappers. The image
has not been built yet; the folder has been checked with oncoanalyser's stub run.

## Contents

| Command | What it is |
|---|---|
| `sage-rs`, `redux-rs`, `esvee-rs` | the ports (`sage-rs/`, `redux-rs/`, `esvee-rs/`) |
| `sage`, `redux`, `esvee` | wrappers (`bin/`) that take the Java tools' command lines as oncoanalyser writes them |
| `samtools`, `sambamba` | given to REDUX and ESVEE as `-bamtool`; the ports only check that the path exists (they sort and index BAMs themselves), so on the host these must be on `PATH` too |

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
