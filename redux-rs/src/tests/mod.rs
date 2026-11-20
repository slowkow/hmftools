//! Ports of the Java unit tests of REDUX 2.1 (hmftools redux/src/test/java/com/hartwig/hmftools/redux).
//! Each test names the Java test it ports. Tests of features the port does not have are present as ignored
//! stubs that give the reason.

pub mod util;

mod bam_merge_test;
mod consensus_read_utils_test;
mod consensus_reads_test;
mod duplicate_group_test;
mod fragment_utils_test;
mod indel_consensus_reads_test;
mod partition_reader_test;
mod read_cache_test;
mod redux_application_test;
mod sorted_bam_writer_test;
mod unmap_reads_test;
