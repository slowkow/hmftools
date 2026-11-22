//! AssemblyStatsTest
use super::*;
use crate::assembly::junction_read_types::calc_proximate_junction_read_ratio;

#[test]
fn test_assembly_supplementary_stats() {
    let _env = test_env();
    let neg_junction = new_junction(CHR_1, 100, REVERSE);
    let assembly_bases = &REF_BASES_200[0..100];
    let base_quals = build_default_base_quals(assembly_bases.len());
    let assembly = JunctionAssembly::new_test(neg_junction, assembly_bases.as_bytes(), &base_quals, 50);
    let read_start_pos = 100;
    // the third supplementary is below the minimum alignment score
    for supp in ["1,1000,+,40M60S,30,0", "1,1400,+,40M60S,30,0", "1,900,+,40M60S,25,0", "2,1000,+,40M60S,30,0"] {
        let read = create_read(TEST_READ_ID, read_start_pos, assembly_bases, "40S60M");
        read.m().bam_record_mut().set_string_attribute(SUPPLEMENTARY_ATTRIBUTE, supp);
        assembly.m().add_junction_read(&read);
    }
    for support_read in assembly.support().clone() {
        assembly.m().stats_mut().add_read(&support_read, &neg_junction, support_read.cached_read().map(|r| &**r));
    }
    assert_close(0.67, assembly.stats().supp_remote_region_ratio(), 0.1);
}

#[test]
fn test_assembly_prox() {
    let _env = test_env();
    let neg_junction = new_junction(CHR_1, 100, REVERSE);
    let assembly_bases = &REF_BASES_200[0..100];
    let raw_reads = vec![
        create_read(TEST_READ_ID, 100, assembly_bases, "40S60M"), // exact and long enough
        create_read(TEST_READ_ID, 100, assembly_bases, "35S65M"), // as above
        create_read(TEST_READ_ID, 100, assembly_bases, "20S80M"), // too short
        create_read(TEST_READ_ID, 70, assembly_bases, "10S90M"),  // too far away to be counted
        create_read(TEST_READ_ID, 90, assembly_bases, "40S60M"),  // included but not exact
    ];
    assert_close(0.5, calc_proximate_junction_read_ratio(&neg_junction, &raw_reads), 0.1);
}
