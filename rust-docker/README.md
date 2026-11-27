# Running the Rust ports in oncoanalyser

The Rust ports of SAGE, REDUX, ESVEE, AMBER, COBALT, the BamTools metrics, TEAL, LILAC,
CIDER, PAVE, PURPLE, LINX, SIGS, Virus Interpreter and PEACH, packaged for use in
place of the Java tools in oncoanalyser
(https://github.com/slowkow/oncoanalyser, branch `rust-ports`, profile `rust_ports`). Two forms:

**A folder on the host**, for running these steps without a container (oncoanalyser's
`--rust_ports_bin`). The ports need nothing beyond the operating system. This builds them and copies
them with the wrappers into one folder:

```bash
rust-docker/install.sh ~/rust-ports
```

**A container image**, for executors that run every task in a container (Google Cloud Batch, AWS
Batch, Kubernetes). Batch machines are x86-64, so build it for linux/amd64. On Google Cloud, Cloud
Build does that and pushes the image to Artifact Registry, without Docker on your machine:

```bash
rust-docker/cloud-build.sh <region>-docker.pkg.dev/<project>/<repository>/hmftools-rust
```

It uploads only the committed files the image needs (about 6 MB; no build folders, nothing
uncommitted), tags the image with the commit id, and takes about 3 minutes on the 32-CPU builder set
in `cloudbuild.yaml`. With Docker, from the repository root:

```bash
docker build --platform linux/amd64 -f rust-docker/Dockerfile -t hmftools-rust:dev .
```

`Dockerfile.dockerignore` limits the build context to the crates, the Java test resources their tests
read, and the wrappers. The build
runs every crate's unit tests on Linux and fails if one fails; the first Linux build found a crash
that macOS hid (rust-htslib's `fetch_seq` freeing htslib's buffer with the Rust allocator), now fixed.
For AWS, build the same Dockerfile and push it to ECR.

## Contents

| Command | What it is |
|---|---|
| `sage-rs`, `redux-rs`, `esvee-rs`, `amber-rs`, `cobalt-rs`, `bam-tools-rs`, `teal-rs`, `lilac-rs`, `cider-rs`, `pave-rs`, `purple-rs`, `linx-rs`, `sigs-rs`, `virus-interpreter-rs`, `peach-rs` | the ports (the folders of the same names) |
| `sage`, `redux`, `esvee`, `amber`, `cobalt`, `bamtools`, `teal`, `lilac`, `cider`, `pave`, `purple`, `linx`, `sigs`, `virusinterpreter`, `peach` | wrappers (`bin/`) that take the Java tools' command lines as oncoanalyser writes them |
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
| `com.hartwig.hmftools.amber.AmberApplication` (or none) | `amber-rs` |
| `com.hartwig.hmftools.cobalt.CobaltApplication` (or none) | `cobalt-rs` |
| `com.hartwig.hmftools.bamtools.metrics.BamMetrics` (or none) | `bam-tools-rs` |
| `com.hartwig.hmftools.teal.TealPipelineTelbamApp`, `TealPipelineApp`, `TealApplication` | `teal-rs` with the class (teal-rs takes the classes itself) |
| `com.hartwig.hmftools.lilac.LilacApplication` (or none) | `lilac-rs` |
| `com.hartwig.hmftools.cider.CiderApplication` (or none) | `cider-rs` (cider-rs takes the class itself) |
| `com.hartwig.hmftools.pave.PaveApplication` (or none) | `pave-rs` |
| `com.hartwig.hmftools.purple.PurpleApplication` (or none) | `purple-rs` (a `-circos` without a value, as `$(which circos)` gives where Circos is missing, is dropped) |
| `com.hartwig.hmftools.linx.LinxApplication` (or none) | `linx-rs` (not the visualiser, `linx.visualiser.SvVisualiser`) |
| `com.hartwig.hmftools.sigs.fitter.SignatureFitter` (or none) | `sigs-rs` |
| `com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication` (or none) | `virus-interpreter-rs` |
| `com.hartwig.hmftools.peach.PeachApplication` (or none) | `peach-rs` |

Any other main class is refused with exit code 2. `-version` prints the version in Java's form, for
example `Sage version 5.1`.

`-java_compat` (or `<TOOL>_JAVA_COMPAT=1`, for example `SAGE_JAVA_COMPAT=1`) makes each port reproduce
the Java tool's output exactly; by default the ports correct the faults listed in their
REPORT.md files.
