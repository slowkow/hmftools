# hmf-common-rs

Code shared by the Rust ports of hmftools (`sage-rs`, `redux-rs`, `esvee-rs`, and the ports that follow).
It plays the part of the Java `hmf-common` module, plus the Java library behaviour the ports reproduce
so that their output equals the Java tools'.

| Module | What it is |
|---|---|
| `args` | command lines in the style of the Java ConfigBuilder (`-name value`, `-flag`) |
| `compat` | the switch between a port's corrected behaviour and Java's (`-java_compat`, `-java_compat_except`, `<TOOL>_JAVA_COMPAT`); each port keeps its own finding ids |
| `faidx`, `refgenome` | indexed FASTA access, 1-based inclusive positions as the Java RefGenomeInterface |
| `fmt` | numbers as Java's `String.format` and `Double.toString` write them |
| `log` | log lines in the layout of the Java tools; the `info!`, `warn!`, `error!` and `debug!` macros |
| `regions` | chromosome ranks and names, regions, `-specific_regions` / `-specific_chr` |
| `java_hash`, `java_id_hash` | java.util.HashMap iteration order, for hash codes and for identity hashes under `-XX:hashCode=2` |
| `java_sort`, `java_math` | TimSort as `List.sort`, and Java numeric behaviour such as `DoubleStream.sum` |

The ports take it as a path dependency (`hmf-common = { path = "../hmf-common-rs" }`), so each keeps
its own `Cargo.lock` and build folder. `cargo test --release` here runs the shared tests; each port's
regression suite is the check for a change here.
