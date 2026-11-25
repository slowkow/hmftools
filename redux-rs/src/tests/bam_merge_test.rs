//! BamMergeTest.

/// BamMergeTest.testSequenceInfos
#[test]
#[ignore = "not ported: tests hmf-common BamMerger.formSequenceIntervals (the merge of the sorted interim BAMs); the port \
            joins its per-partition block streams directly and has no counterpart"]
fn test_sequence_infos() {}

// not in Java: finding R41, the order of the input BAMs

/// A header-only BAM with one read group, in the test folder
fn lane_bam(name: &str, read_group: &str) -> String {
    use rust_htslib::bam::{header::HeaderRecord, Format, Header, Writer};
    let dir = std::env::temp_dir().join(format!("redux_rs_r41_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name).to_string_lossy().to_string();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "1").push_tag(b"LN", 1000));
    header.push_record(HeaderRecord::new(b"RG").push_tag(b"ID", read_group).push_tag(b"SM", "TEST"));
    drop(Writer::from_path(&path, &header, Format::Bam).unwrap());
    path
}

fn config_inputs(inputs: &[&String]) -> Vec<String> {
    let joined: Vec<&str> = inputs.iter().map(|s| s.as_str()).collect();
    let args: Vec<String> = ["-sample", "TEST", "-input_bam", &joined.join(","), "-ref_genome", super::util::placeholder_fasta()]
        .iter()
        .map(|s| s.to_string())
        .collect();
    crate::config::Config::from_args(&args).unwrap().bam_files
}

#[test]
fn corrected_inputs_are_taken_in_header_order_whatever_the_order_given() {
    use crate::compat::with_java;
    // the file names sort the other way round from the read groups
    let a = lane_bam("a.lane.bam", "LANE_2");
    let b = lane_bam("b.lane.bam", "LANE_1");
    let c = lane_bam("c.lane.bam", "LANE_3");
    with_java(false, || {
        let expected = vec![b.clone(), a.clone(), c.clone()];
        for given in [[&a, &b, &c], [&c, &b, &a], [&b, &c, &a]] {
            assert_eq!(expected, config_inputs(&given));
        }
        // one input is left as it is
        assert_eq!(vec![a.clone()], config_inputs(&[&a]));
    });
    // Java: the order given
    with_java(true, || {
        assert_eq!(vec![a.clone(), b.clone(), c.clone()], config_inputs(&[&a, &b, &c]));
        assert_eq!(vec![c.clone(), b.clone(), a.clone()], config_inputs(&[&c, &b, &a]));
    });
}

#[test]
fn inputs_with_the_same_header_keep_the_order_of_their_paths() {
    let a = lane_bam("same_a.bam", "LANE");
    let b = lane_bam("same_b.bam", "LANE");
    let ordered = crate::config::input_order_by_header(vec![b.clone(), a.clone()]).unwrap();
    assert_eq!(vec![a, b], ordered);
}
