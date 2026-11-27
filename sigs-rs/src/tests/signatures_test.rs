//! SignaturesTest (sigs/src/test/java/com/hartwig/hmftools/sigs/SignaturesTest.java)

use hmf_common::variant_context::{VcfHeader, VcfRecord};

use crate::config::VariantFilters;
use crate::snv_loader::{create_variant, SigVariant};

const TEST_SAMPLE_ID: &str = "Sample_001";

fn header() -> VcfHeader {
    VcfHeader::from_lines(
        vec![
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">".into(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">".into(),
            "##INFO=<ID=HOTSPOT,Number=0,Type=Flag,Description=\"Hotspot\">".into(),
        ],
        vec![TEST_SAMPLE_ID.to_string()],
    )
}

/// buildVariantContext(variant, altCount, depth) with extra INFO fields
fn build(alt_count: i32, depth: i32, info: &str) -> SigVariant {
    let h = header();
    let line = format!("chr1\t100\t.\tA\tT\t.\t.\t{info}\tGT:AD:DP\t./.:0,{alt_count}:{depth}");
    let mut rec = VcfRecord::parse(&line, &h).unwrap();
    create_variant(&mut rec, &h, TEST_SAMPLE_ID).unwrap().unwrap()
}

#[test]
fn test_variant_filters() {
    let filters = VariantFilters { tumor_only: true, ..Default::default() };

    // AF 4% is below the tumor-only threshold
    assert!(!filters.passes(&build(4, 100, ".")));

    let variant = build(10, 100, ".");
    assert!(filters.passes(&variant));

    // a non-zero gnomAD frequency
    let mut with_gnomad = variant.clone();
    with_gnomad.gnomad_freq = 0.001;
    assert!(!filters.passes(&with_gnomad));
    assert!(filters.passes(&variant));

    assert!(!filters.passes(&build(10, 100, "HOTSPOT")));
    assert!(!filters.passes(&build(10, 100, "SOM_LH=LOW")));
}
