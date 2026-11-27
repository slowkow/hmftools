//! With the `bwa` feature: builds the vendored BWA-MEM alignment library (vendor/bwa, see its README.md)
//! for the `bwa` module. Without it, nothing.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    #[cfg(feature = "bwa")]
    {
        let dir = "vendor/bwa";
        let sources = [
            "utils.c", "kthread.c", "kstring.c", "ksw.c", "bwt.c", "bntseq.c", "bwa.c", "bwamem.c", "bwamem_pair.c", "bwamem_extra.c", "malloc_wrap.c",
            "hmf_bwa.c",
        ];
        let mut build = cc::Build::new();
        for s in sources {
            build.file(format!("{dir}/{s}"));
            println!("cargo:rerun-if-changed={dir}/{s}");
        }
        build
            .include(dir)
            .define("HAVE_PTHREAD", None)
            .define("USE_MALLOC_WRAPPERS", None)
            // htslib has a function of the same name
            .define("ksprintf", "hmf_bwa_ksprintf")
            .opt_level(2)
            .warnings(false)
            .flag_if_supported("-std=gnu99");
        build.compile("hmf_bwa");
        println!("cargo:rustc-link-lib=z");
    }
}
