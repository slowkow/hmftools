//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/candidate/VariantTypesTest.java
//!
//! The Java tests run a RegionTask through RegionTaskTester, with a MockRefGenome and a MockSamSlicer
//! holding the tumor reads. The port's region task (`process_region` then `finalise_region` in main.rs)
//! reads an indexed FASTA and an indexed BAM, so `RegionTaskTester` here writes both as temporary files.
//!
//! RegionTask.getVariants() is every variant that passed the hard filters, whatever its soft filters;
//! that is `RegionOutput.variants` as returned by `process_region`. `finalise_region` (phasing, dedup,
//! final filters) is then run as well, as RegionTask.run() does, but nothing is asserted on its output.

use super::candidate_creation::generate_random_bases;
use super::indel_dedup::{mock_ref_genome, test_tmp_dir};
use super::util::*;
use crate::bamio::BamSlicer;
use crate::quality::{BqrMap, MsiJitter};
use crate::read_context::VariantReadContext;
use crate::types::*;
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use std::sync::Arc;

// "X" covers the zero-position
// pos:    1         11        21        31        41       50 and repeats up to 150
const TEST_REF_BASES: &str = concat!(
    "X",
    "GCAGGAGAATCCCTTGAACCTGGGGGGAACCTGGGGGGGTGAGCTGAGAT",
    "CATGCCATTGCACTCTAGCCTGGGCAACAAGAGTGAAACTCCGCCTCAAA",
    "ACAAACAAACAAACAAACAAACAAACAAACAAACAAAAACCTCCAAAACA"
);

/// What the Java tests read from a SageVariant
struct TestVariant {
    position: i32,
    ref_: String,
    alt: String,
    read_context: Arc<VariantReadContext>,
}

/// RegionTaskTester: default config, tumor sample only, no hotspots, panel, transcripts or high
/// confidence regions, no BQR records and default jitter parameters.
struct RegionTaskTester {
    test_name: &'static str,
    /// RefGenome.RefGenomeMap.get(CHR_1)
    ref_bases: String,
    /// TumorSamSlicer.ReadRecords
    tumor_reads: Vec<Read>,
}

impl RegionTaskTester {
    fn new(test_name: &'static str) -> Self {
        RegionTaskTester { test_name, ref_bases: String::new(), tumor_reads: Vec::new() }
    }

    fn write_tumor_bam(&self) -> String {
        let path = test_tmp_dir().join(format!("{}.bam", self.test_name));
        let path = path.to_str().unwrap().to_string();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", CHR_1);
        sq.push_tag(b"LN", self.ref_bases.len() - 1);
        header.push_record(&sq);

        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            for read in &self.tumor_reads {
                let cigar = CigarString(
                    read.cigar
                        .iter()
                        .map(|el| {
                            let len = el.len as u32;
                            match el.op {
                                Op::M => Cigar::Match(len),
                                Op::I => Cigar::Ins(len),
                                Op::D => Cigar::Del(len),
                                Op::N => Cigar::RefSkip(len),
                                Op::S => Cigar::SoftClip(len),
                                Op::H => Cigar::HardClip(len),
                                Op::P => Cigar::Pad(len),
                                Op::Eq => Cigar::Equal(len),
                                Op::X => Cigar::Diff(len),
                            }
                        })
                        .collect(),
                );
                let mut record = bam::Record::new();
                record.set(read.name.as_bytes(), Some(&cigar), &read.bases, &read.quals);
                record.set_tid(read.tid);
                record.set_pos(read.start as i64 - 1);
                record.set_mapq(read.mapq as u8);
                record.set_flags(read.flags);
                record.set_mtid(read.mate_tid);
                record.set_mpos(read.mate_start as i64 - 1);
                record.set_insert_size(read.insert_size as i64);
                if let Some(nm) = read.nm {
                    record.push_aux(b"NM", Aux::I32(nm)).unwrap();
                }
                if let Some(mate_cigar) = &read.mate_cigar {
                    record.push_aux(b"MC", Aux::String(mate_cigar)).unwrap();
                }
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        path
    }

    /// createRegionTask(region), task.run() and task.getVariants()
    fn run_region_task(&self, region: BaseRegion) -> Vec<TestVariant> {
        let cfg = test_config();
        let ref_data = empty_ref_data();
        let bqr = [BqrMap::default()];
        let msi = [MsiJitter { params: Vec::new(), probable_msi: false }];
        let shared = crate::Shared { cfg: &cfg, ref_data: &ref_data, bqr: &bqr, msi: &msi };

        let genome = mock_ref_genome(self.test_name, &[(CHR_1, &self.ref_bases)]);
        let tumor_bam = BamSlicer::open(&self.write_tumor_bam(), 1).unwrap();
        let mut worker = crate::Worker { genome, tumor_bam, ref_bam: None };

        let Some(output) = crate::process_region(&shared, &mut worker, CHR_1, region).unwrap() else {
            return Vec::new();
        };

        let variants: Vec<TestVariant> = output
            .variants
            .iter()
            .map(|v| TestVariant {
                position: v.position(),
                ref_: v.variant().ref_str().to_string(),
                alt: v.variant().alt_str().to_string(),
                read_context: v.candidate.ctx.clone(),
            })
            .collect();

        if !variants.is_empty() {
            let mut next_lps = 1;
            crate::finalise_region(&shared, &worker.genome, CHR_1, region, output, &mut next_lps);
        }

        variants
    }
}

/// The Java test gives both of its reads 51 bases with a 50M cigar. MockSamSlicer hands such a record
/// straight to SAGE, but htslib refuses to read one back from a BAM, so the port cannot be given it. The
/// 51st base lies beyond the cigar and is never aligned; it is dropped here so the record is valid. This
/// is the only change to the Java inputs.
fn valid_bam_read(mut read: Read) -> Read {
    assert_eq!(51, read.bases.len());
    let cigar_read_length: i32 = read.cigar.iter().filter(|el| el.op.consumes_read()).map(|el| el.len).sum();
    read.bases.truncate(cigar_read_length as usize);
    read.quals.truncate(cigar_read_length as usize);
    read
}

// VariantTypesTest.testBasicVariants
#[test]
fn test_basic_variants() {
    let region = BaseRegion::new(1, 200);

    let mut tester = RegionTaskTester::new("variant_types_basic_variants");

    let ref_bases = TEST_REF_BASES;
    tester.ref_bases = format!("{}{}", ref_bases, generate_random_bases(1500)); // need to cover the ref sequence buffer

    let mut read_bases = format!("{}A{}", &ref_bases[1..21], &ref_bases[22..52]);

    let read1 = create_sam_record("READ_01", 1, &read_bases, "50M");
    tester.tumor_reads.push(valid_bam_read(read1));

    read_bases = format!("{}A{}", &ref_bases[3..21], &ref_bases[22..54]);
    let read2 = create_sam_record("READ_02", 3, &read_bases, "50M");
    tester.tumor_reads.push(valid_bam_read(read2));

    let variants = tester.run_region_task(region);

    assert_eq!(1, variants.len());
    let var = &variants[0];
    assert_eq!(21, var.position);
    assert_eq!("T", var.ref_);
    assert_eq!("A", var.alt);
}

// VariantTypesTest.testInsertWithRepeats
#[test]
fn test_insert_with_repeats() {
    let region = BaseRegion::new(1, 200);

    let mut tester = RegionTaskTester::new("variant_types_insert_with_repeats");

    let ref_bases = TEST_REF_BASES;
    tester.ref_bases = format!("{}{}", ref_bases, generate_random_bases(1500)); // need to cover the ref sequence buffer

    let inserted_bases = "GGGGGGGG";
    let read_start_pos = 81;

    // extend the read to the end of the sequence so it's past the section of repeats
    let read_bases = format!("{}{}{}", &ref_bases[read_start_pos..111], inserted_bases, &ref_bases[111..151]);

    let read1 = create_sam_record("READ_01", read_start_pos as i32, &read_bases, "30M8I40M");
    let read2 = create_sam_record("READ_02", read_start_pos as i32, &read_bases, "30M8I40M");
    let read3 = create_sam_record("READ_03", read_start_pos as i32, &read_bases, "30M8I40M");
    tester.tumor_reads.push(read1);
    tester.tumor_reads.push(read2);
    tester.tumor_reads.push(read3);

    let variants = tester.run_region_task(region);

    assert_eq!(1, variants.len());
    let var = &variants[0];
    assert_eq!(110, var.position);
    assert_eq!("C", var.ref_);
    assert_eq!(format!("C{inserted_bases}"), var.alt);

    let read_context = &var.read_context;
    assert!(read_context.is_valid());
    assert_eq!(49, read_context.core_length());
    assert!(read_context.max_repeat.is_some());
    let max_repeat = read_context.max_repeat.as_ref().unwrap();
    assert_eq!("G", max_repeat.bases_str());
    assert_eq!(8, max_repeat.count);
}

// VariantTypesTest.testInsertWithHomology
#[test]
fn test_insert_with_homology() {
    let region = BaseRegion::new(1, 200);

    let mut tester = RegionTaskTester::new("variant_types_insert_with_homology");

    let ref_bases = TEST_REF_BASES;
    tester.ref_bases = format!("{}{}", ref_bases, generate_random_bases(1500)); // need to cover the ref sequence buffer

    let inserted_bases = "GGGGGGAACCT";
    let read_start_pos = 1;
    let var_position = 34;

    // extend the read to the end of the sequence so it's past the section of repeats
    let read_bases =
        format!("{}{}{}", &ref_bases[read_start_pos..var_position + 1], inserted_bases, &ref_bases[var_position + 1..81]);

    let read1 = create_sam_record("READ_01", read_start_pos as i32, &read_bases, "34M11I46M");
    let read2 = create_sam_record("READ_02", read_start_pos as i32, &read_bases, "34M11I46M");
    let read3 = create_sam_record("READ_03", read_start_pos as i32, &read_bases, "34M11I46M");
    tester.tumor_reads.push(read1);
    tester.tumor_reads.push(read2);
    tester.tumor_reads.push(read3);

    let variants = tester.run_region_task(region);

    assert_eq!(1, variants.len());
    let var = &variants[0];

    // no left alignment
    assert_eq!(var_position as i32, var.position);
    assert_eq!("G", var.ref_);
    assert_eq!(format!("G{inserted_bases}"), var.alt);
    let read_context = &var.read_context;
    assert!(read_context.is_valid());
    assert_eq!(21, read_context.core_length());
    let homology_bases = read_context.homology.as_ref().map(|h| String::from_utf8(h.bases.clone()).unwrap());
    assert_eq!(Some("GGGGG".to_string()), homology_bases);
}
