//! Builds the vendored BWA-MEM alignment library used by the assembly step (see vendor/bwa/README.md).

fn main() {
    let dir = "vendor/bwa";
    let sources = [
        "utils.c", "kthread.c", "kstring.c", "ksw.c", "bwt.c", "bntseq.c", "bwa.c", "bwamem.c", "bwamem_pair.c", "bwamem_extra.c", "malloc_wrap.c", "esvee_bwa.c",
    ];
    let mut build = cc::Build::new();
    for s in sources {
        build.file(format!("{dir}/{s}"));
        println!("cargo:rerun-if-changed={dir}/{s}");
    }
    build.include(dir).define("HAVE_PTHREAD", None).define("USE_MALLOC_WRAPPERS", None)
        // htslib has a function of the same name
        .define("ksprintf", "esvee_bwa_ksprintf").opt_level(2).warnings(false).flag_if_supported("-std=gnu99");
    build.compile("esvee_bwa");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rerun-if-changed=build.rs");
}
