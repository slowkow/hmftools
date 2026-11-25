//! Not from Java: the corrected behaviour for the findings of REPORT.md section 5, beside the Java
//! behaviour that `-java_compat` keeps. Phasing (S1 to S3) is tested in `phase`, TINC (S10 to S12) in
//! `tinc`, the append records (S19, S20) in `append`, the CRAM NM tag (S7) in `alignment_reading`.

use super::util::*;
use crate::compat::with_java;
use crate::types::*;
use std::io::Write;

fn temp_file(name: &str, content: &str) -> String {
    let path = std::env::temp_dir().join(format!("sage_rs_corrections_{}_{name}", std::process::id()));
    std::fs::File::create(&path).unwrap().write_all(content.as_bytes()).unwrap();
    path.to_string_lossy().to_string()
}

fn config_with(extra: &[&str]) -> crate::config::Config {
    let mut args: Vec<String> =
        ["-tumor", "TUMOR_ID", "-ref_genome", "none.fasta", "-ref_genome_version", "37", "-output_vcf", "none.vcf.gz"]
            .iter()
            .map(|x| x.to_string())
            .collect();
    args.extend(extra.iter().map(|x| x.to_string()));
    crate::config::Config::from_args(&args).unwrap()
}

/// S4: the support is scaled by the quality ratio, then rounded
#[test]
fn vaf_probability_support_is_scaled_before_rounding() {
    use crate::filters::adjusted_alt_support_count;
    assert_eq!(adjusted_alt_support_count(10, 0.74), 7);
    assert_eq!(adjusted_alt_support_count(10, 0.4), 4);
    assert_eq!(adjusted_alt_support_count(10, 1.3), 10);
    assert_eq!(with_java(true, || adjusted_alt_support_count(10, 0.74)), 10);
    assert_eq!(with_java(true, || adjusted_alt_support_count(10, 0.4)), 0);
}

/// S5: a pair whose merge fails is closed; the first mate is processed at once and nothing stays cached
#[test]
fn failed_fragment_merge_leaves_nothing_cached() {
    let run = || {
        let bases = "ACGTACGTACGTACGTACGT";
        let mut first = create_sam_record("READ_X", 100, bases, "20M");
        first.mate_start = 105;
        let mut second = create_sam_record("READ_X", 105, bases, "20M");
        second.flags = FLAG_PAIRED | FLAG_PROPER | FLAG_REVERSE;
        second.mate_start = 100;
        second.cigar.clear(); // the merge cannot walk an empty cigar: Java throws
        let mut sync = crate::sync::FragmentSync::new();
        let pool = crate::bamio::ReadPool::default();
        let mut out = Vec::new();
        assert!(sync.handle(first, &mut out, &pool).is_none());
        assert!(out.is_empty());
        assert!(sync.handle(second, &mut out, &pool).is_none());
        assert_eq!(sync.counts[crate::sync::SYNC_EXCEPTION], 1);
        let processed: Vec<i32> = out.iter().map(|(read, _)| read.start).collect();
        (processed, sync.empty_cached_reads().len())
    };
    assert_eq!(run(), (vec![100], 0));
    assert_eq!(with_java(true, run), (vec![], 1));
}

/// S6: a BED file that comes back to a chromosome keeps both blocks, in position order
#[test]
fn bed_file_keeps_every_block_of_a_chromosome() {
    let bed = temp_file("revisit.bed", "1\t500\t600\n2\t10\t20\n1\t100\t200\n1\t700\t800\n");
    let regions = |map: &std::collections::HashMap<String, Vec<BaseRegion>>, chr: &str| -> Vec<(i32, i32)> {
        map[chr].iter().map(|r| (r.start, r.end)).collect()
    };
    let map = crate::refdata::load_bed_chr_map(&bed, false).unwrap();
    assert_eq!(regions(&map, "1"), vec![(101, 200), (501, 600), (701, 800)]);
    assert_eq!(regions(&map, "2"), vec![(11, 20)]);
    let java = with_java(true, || crate::refdata::load_bed_chr_map(&bed, false).unwrap());
    assert_eq!(regions(&java, "1"), vec![(101, 200), (701, 800)]);
    let _ = std::fs::remove_file(bed);
}

/// S15: the Ultima maximum raw quality is the maximum over the samples
#[test]
fn ultima_max_raw_qual_is_the_maximum_over_samples() {
    let table = |qual: i32| {
        temp_file(
            &format!("bqr_{qual}.tsv"),
            &format!("Alt\tRef\tTrinucleotideContext\tOriginalQual\tRecalibratedQual\nA\tC\tACG\t{qual}\t30.0\nA\tC\tACG\t20\t18.0\n"),
        )
    };
    let files = [table(35), table(40), table(37)];
    let bqr: Vec<crate::quality::BqrMap> = files.iter().map(|f| crate::quality::BqrMap::load(f).unwrap()).collect();
    for (tumors, references) in [("T1", "N1,N2"), ("SAMPLE_B", "SAMPLE_A,SAMPLE_C"), ("X", "Z,Y")] {
        let cfg = config_with(&["-tumor", tumors, "-reference", references, "-reference_bam", "a.bam,b.bam", "-tumor_bam", "t.bam"]);
        assert_eq!(crate::ultima_max_raw_qual(&cfg, &bqr), 40);
        // Java: the value of whichever sample comes last in its map
        assert!([35, 37, 40].contains(&with_java(true, || crate::ultima_max_raw_qual(&cfg, &bqr))));
    }
    let cfg = config_with(&["-tumor", "T1", "-tumor_bam", "t.bam", "-skip_bqr"]);
    assert_eq!(crate::ultima_max_raw_qual(&cfg, &bqr[..1]), 0);
    for f in files {
        let _ = std::fs::remove_file(f);
    }
}

/// S24: a short contig is cut the same way at every thread count
#[test]
fn partition_does_not_depend_on_threads() {
    let regions = |threads: &str| -> Vec<(i32, i32)> {
        let cfg = config_with(&["-threads", threads]);
        crate::partition(&cfg, "MT", 16569).iter().map(|r| (r.start, r.end)).collect()
    };
    assert_eq!(regions("1"), vec![(1, 16569)]);
    assert_eq!(regions("8"), regions("1"));
    assert_eq!(with_java(true, || regions("1")), vec![(1, 16569)]);
    assert_eq!(with_java(true, || regions("8")).len(), 4);
    // a long chromosome is cut into slices of the configured size in both modes
    let cfg = config_with(&["-threads", "8"]);
    let long = crate::partition(&cfg, "1", 1_000_000);
    assert_eq!(long.len(), 10);
    assert_eq!(with_java(true, || crate::partition(&cfg, "1", 1_000_000)).len(), 10);
}

/// S38: `Ter` is the stop in a substitution; an unknown amino acid code is an error
#[test]
fn protein_change_codes() {
    use crate::vis_aa::{parse_amino_acid_events, AminoAcidEvent};
    assert_eq!(parse_amino_acid_events("p.Gly12Ter").unwrap(), vec![AminoAcidEvent::Substitution { pos: 12, alt: b'*' }]);
    assert_eq!(parse_amino_acid_events("p.Gly12*").unwrap(), vec![AminoAcidEvent::Substitution { pos: 12, alt: b'*' }]);
    assert!(with_java(true, || parse_amino_acid_events("p.Gly12Ter")).is_err());
    assert!(parse_amino_acid_events("p.Gly12_Ala13delinsXyzVal").is_err());
    // Java: the letters of "null" and V, five amino acids for a gap of two
    let java = with_java(true, || parse_amino_acid_events("p.Gly12_Ala13delinsXyzVal")).unwrap();
    assert_eq!(java.len(), 3);
    assert!(parse_amino_acid_events("p.Gly12Xyz").is_err());
    assert!(parse_amino_acid_events("p.Gly12insAla").is_err());
}

/// S40, S41: the header declares what is written
#[test]
fn header_declarations() {
    let rc_jit = |lines: Vec<&'static str>| lines.into_iter().find(|l| l.starts_with("##FORMAT=<ID=RC_JIT,")).unwrap();
    assert!(rc_jit(crate::vcfout::header_field_lines().collect()).contains("Number=2,"));
    assert!(with_java(true, || rc_jit(crate::vcfout::header_field_lines().collect())).contains("Number=3,"));

    let header = |name: &str| -> String {
        let path = std::env::temp_dir().join(format!("sage_rs_corrections_{}_{name}.vcf.gz", std::process::id()));
        let path = path.to_string_lossy().to_string();
        let writer =
            crate::vcfout::VcfWriter::create(&path, &[("1".to_string(), 1000)], &[], &["TUMOR_ID".to_string()], false).unwrap();
        writer.close().unwrap();
        let mut text = String::new();
        use std::io::Read as _;
        flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap()).read_to_string(&mut text).unwrap();
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(format!("{path}.tbi"));
        text
    };
    let corrected = header("corrected");
    assert!(!corrected.contains("##FILTER=<ID=fragmentLength,"));
    assert!(corrected.contains("##FILTER=<ID=minTumorQual,"));
    let java = with_java(true, || header("java"));
    assert!(java.contains("##FILTER=<ID=fragmentLength,"));
    assert_eq!(java.lines().count(), corrected.lines().count() + 1);
}

/// S44: the overlap of the mates is marked on the columns they share, also after an insert
#[test]
fn overlap_of_mates_after_an_insert() {
    use crate::vis_svg::BaseSeqViewModel;
    let merged = build_sam_record(100, "4M2I4M", "ACGTTTACGT");
    let first = BaseSeqViewModel::from_read(&build_sam_record(104, "4M", "ACGT"));
    let second = BaseSeqViewModel::from_read(&build_sam_record(100, "8M", "ACGTACGT"));
    let overlapped = || -> Vec<bool> {
        BaseSeqViewModel::from_consensus_fragment(&merged, Some(&first), Some(&second)).bases.iter().map(|b| b.overlapped).collect()
    };
    assert_eq!(overlapped(), vec![false, false, false, false, true, true, true, true]);
    assert_eq!(with_java(true, overlapped), vec![false, false, false, false, true, true, false, false]);
}

/// S49: a deletion of more than 5 repeat units is beyond the jitter model, as an insertion is
#[test]
fn msi_error_rate_is_zero_beyond_five_units_either_way() {
    use crate::quality::{JitterParams, MsiJitter};
    use crate::repeat::RepeatInfo;
    let params = JitterParams {
        repeat_unit: "A/T".to_string(), consensus: "NONE".to_string(), scale4: 0.05, scale5: 0.06, scale6: 0.07, gradient: 0.0444,
        intercept: -0.1835, skew: 1.0164,
    };
    let msi = MsiJitter { params: vec![params], probable_msi: false };
    let repeat = RepeatInfo::new(1, b"A".to_vec(), 12);
    let del6 = SimpleVariant::new(100, b"GAAAAAA", b"G");
    let ins6 = SimpleVariant::new(100, b"G", b"GAAAAAA");
    let del5 = SimpleVariant::new(100, b"GAAAAA", b"G");
    with_java(false, || {
        assert_eq!(0.0, msi.calc_error_rate(&del6, &repeat));
        assert_eq!(0.0, msi.calc_error_rate(&ins6, &repeat));
        assert!(msi.calc_error_rate(&del5, &repeat) > 0.0);
    });
    // Java: the 6-unit deletion gets the model's floor, the 6-unit insertion none
    with_java(true, || {
        assert!(msi.calc_error_rate(&del6, &repeat) > 0.0);
        assert_eq!(0.0, msi.calc_error_rate(&ins6, &repeat));
    });
}

/// S50: an MNV filtered by an earlier pair does not filter another one
#[test]
fn filtered_mnv_does_not_filter_the_next() {
    use super::variant_dedup::{add_local_phase_set, create_read_context, create_sage_variant_from_context, set_tumor_quality};
    use crate::dedup::dedup_mnv_overlaps;
    use crate::filters::SoftFilter;
    let mnv = |pos: i32, r: &str, a: &str, count: i32, qual: i32| {
        let v = SimpleVariant::new(pos, r.as_bytes(), a.as_bytes());
        let mut s = create_sage_variant_from_context(create_read_context(v, "A", "A", "GGGGGG", "GGGGGG"));
        set_tumor_quality(&mut s, count, qual);
        add_local_phase_set(&mut s, 1, count);
        s
    };
    // A covers 100-102, B 102-104 and C 103-104: A filters B (equal changes, higher quality); C overlaps only B
    let make = || vec![mnv(100, "ACG", "CAT", 5, 1000), mnv(102, "GTA", "TGC", 4, 900), mnv(103, "TA", "GC", 3, 800)];
    let mut corrected = make();
    with_java(false, || dedup_mnv_overlaps(&mut corrected));
    assert!(corrected[0].is_passing());
    assert!(corrected[1].filters.contains(&SoftFilter::DedupMnv));
    assert!(corrected[2].is_passing());
    // Java: the filtered B still filters C (3 changed bases against 2)
    let mut java = make();
    with_java(true, || dedup_mnv_overlaps(&mut java));
    assert!(java[0].is_passing());
    assert!(!java[1].is_passing() && !java[2].is_passing());
}
