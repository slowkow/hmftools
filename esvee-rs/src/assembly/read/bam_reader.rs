//! BamReader: slices the prep BAMs of all samples for a region, tumor files first.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as _};

use crate::assembly::config::AssemblyConfig;
use crate::assembly::read::read_heap;
use crate::assembly::sam::*;

struct SampleReader {
    reader: bam::IndexedReader,
    /// interned contig names by target id
    names: Vec<&'static str>,
}

pub struct BamReader {
    readers: Vec<SampleReader>,
    tumor_count: usize,
    current_is_reference_sample: bool,
    record: bam::Record,
}

const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

impl BamReader {
    pub fn new(config: &AssemblyConfig) -> Result<BamReader> {
        let mut readers = Vec::new();
        for bam_file in config.combined_bam_files() {
            let reader = bam::IndexedReader::from_path(&bam_file).with_context(|| format!("failed to open BAM: {bam_file}"))?;
            let names = reader.header().target_names().iter().map(|n| intern(&String::from_utf8_lossy(n))).collect();
            readers.push(SampleReader { reader, names });
        }
        Ok(BamReader { readers, tumor_count: config.tumor_ids.len(), current_is_reference_sample: false, record: bam::Record::new() })
    }

    pub fn current_is_reference_sample(&self) -> bool {
        self.current_is_reference_sample
    }

    /// Calls the consumer with each record overlapping the region (BamSlicer with unmapped reads and
    /// supplementaries kept; duplicates and secondaries without a hard clip dropped), then whether the
    /// record came from a reference sample.
    pub fn slice_bam<F: FnMut(SamRecord, bool)>(&mut self, config: &AssemblyConfig, chromosome: &str, position_start: i32, position_end: i32, mut consumer: F) {
        let bam_pos_start = position_start.max(1);
        let bam_pos_end = position_end.min(config.chromosome_length(chromosome));
        if bam_pos_start > bam_pos_end {
            return;
        }
        for i in 0..self.readers.len() {
            let is_reference = i >= self.tumor_count;
            self.current_is_reference_sample = is_reference;
            let sr = &mut self.readers[i];
            let Some(tid) = sr.reader.header().tid(chromosome.as_bytes()) else {
                crate::warn!("cannot find sequence index for chromosome {} in bam header", chromosome);
                continue;
            };
            // the Java slice uses the requested coordinates; a start below 1 is the contig start
            if sr.reader.fetch((tid, (position_start.max(1) - 1) as i64, position_end as i64)).is_err() {
                continue;
            }
            while let Some(Ok(())) = sr.reader.read(&mut self.record) {
                let rec = &self.record;
                let flags = rec.flags() as i32;
                if flags & FLAG_SECONDARY != 0 && !rec.cigar().iter().any(|c| matches!(c, bam::record::Cigar::HardClip(_))) {
                    continue;
                }
                if flags & FLAG_DUPLICATE != 0 {
                    continue;
                }
                consumer(to_sam_record(rec, &sr.names, i as i32), is_reference);
            }
        }
    }
}

fn to_sam_record(rec: &bam::Record, names: &[&'static str], sample_index: i32) -> SamRecord {
    let name_of = |tid: i32| if tid >= 0 && (tid as usize) < names.len() { names[tid as usize] } else { NO_CHROMOSOME_NAME };
    // the buffers of the record are allocated in the read heap (see read_heap)
    let raw_cigar = rec.raw_cigar();
    let cigar: Vec<CigarElement> = read_heap::vec_from(&raw_cigar.iter().map(|c| CigarElement::new((c >> 4) as i32, CigarOp::from_bam(c & 0xf))).collect::<Vec<_>>());
    let seq = rec.seq();
    let enc = seq.encoded;
    let len = seq.len();
    let mut read_bases: Vec<u8> = read_heap::vec_from(&rec.qual()[..len]);
    for (i, base) in read_bases.iter_mut().enumerate() {
        let b = enc[i >> 1];
        let code = if i & 1 == 0 { b >> 4 } else { b & 0xf };
        *base = BASES[code as usize];
    }
    let qual = rec.qual();
    let base_qualities = if !qual.is_empty() && qual[0] == 0xff { Vec::new() } else { read_heap::vec_from(qual) };
    // auxiliary data: after the fixed fields, name, cigar, sequence and qualities
    let data = rec.inner().data;
    let l_data = rec.inner().l_data as usize;
    let core = &rec.inner().core;
    let aux_offset = core.l_qname as usize + core.n_cigar as usize * 4 + (len + 1) / 2 + len;
    let aux = if aux_offset < l_data { read_heap::vec_from(&copy_used_attributes(unsafe { std::slice::from_raw_parts(data.add(aux_offset), l_data - aux_offset) })) } else { Vec::new() };
    SamRecord {
        read_name: read_heap::string_from(&String::from_utf8_lossy(rec.qname())),
        flags: rec.flags() as i32,
        reference_name: name_of(rec.tid()),
        alignment_start: rec.pos() as i32 + 1,
        mapping_quality: rec.mapq() as i32,
        cigar,
        mate_reference_name: name_of(rec.mtid()),
        mate_alignment_start: rec.mpos() as i32 + 1,
        inferred_insert_size: rec.insert_size() as i32,
        read_bases,
        base_qualities,
        aux,
        sample_index,
    }
}
