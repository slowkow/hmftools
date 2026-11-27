//! The NEO_EPITOPES function: RNA fragment support for NEO's fusion neo-epitopes (NeoEpitopeReader,
//! NeoEpitopeData, NeoFragmentMatcher, NeoFragmentSupport).

use crate::bam::{to_read, BamReader, SliceFilter};
use crate::common::*;
use crate::config::IsofoxConfig;
use crate::gene::{create_gene_read_data, GeneCollection};
use crate::read::{parse_cigar, Read};
use crate::read_trans::{find_overlapping_regions, process_overlapping_regions};
use crate::refgen::RefSource;
use anyhow::{Context, Result};
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::java_seq_hash::SeqHashMap;
use hmf_common::neo::{read_neo_epitope_file, rna_neo_epitope_header, NeoEpitopeFile};
use hmf_common::regions::BaseRegion;
use rust_htslib::bam;

const EXACT_MATCH: usize = 0;
const PARTIAL_MATCH: usize = 1;
const MISMATCH: usize = 2;
const MIN_BASE_OVERLAP: i32 = 10;
const MAX_BASE_MISMATCH: i32 = 2;

#[derive(Clone, Debug, Default)]
struct Support {
    up: [i32; 2],
    novel: [i32; 2],
    down: [i32; 2],
    ref_base_depth: [i32; 2],
}

impl Support {
    fn combine(&mut self, o: &Support) {
        for i in EXACT_MATCH..=PARTIAL_MATCH {
            self.up[i] += o.up[i];
            self.novel[i] += o.novel[i];
            self.down[i] += o.down[i];
        }
    }
    fn set_max(&mut self, o: &Support) {
        for i in EXACT_MATCH..=PARTIAL_MATCH {
            self.up[i] = self.up[i].max(o.up[i]);
            self.novel[i] = self.novel[i].max(o.novel[i]);
            self.down[i] = self.down[i].max(o.down[i]);
        }
    }
}

/// NeoEpitopeData
struct NeoData {
    source: NeoEpitopeFile,
    chromosomes: [String; 2],
    positions: [i32; 2],
    orientations: [i8; 2],
    transcripts: [Vec<String>; 2],
    coding_base_coords: [Vec<[i32; 2]>; 2],
    support: Support,
}

/// SamRecordUtils.generateMappedCoords
fn generate_mapped_coords(cigar: &[crate::read::CigarElem], pos_start: i32) -> Vec<[i32; 2]> {
    use crate::read::{OP_D, OP_I, OP_M, OP_N, OP_S};
    let mut out: Vec<[i32; 2]> = Vec::new();
    let mut offset = 0;
    let mut cont = false;
    for e in cigar {
        match e.op {
            OP_S => {}
            OP_D => {
                offset += e.len as i32;
                cont = true;
            }
            OP_I => cont = true,
            OP_N => {
                offset += e.len as i32;
                cont = false;
            }
            OP_M => {
                let s = pos_start + offset;
                let end = s + e.len as i32 - 1;
                if cont && !out.is_empty() {
                    out.last_mut().unwrap()[1] = end;
                } else {
                    out.push([s, end]);
                }
                offset += e.len as i32;
                cont = false;
            }
            _ => {}
        }
    }
    out
}

impl NeoData {
    fn new(source: NeoEpitopeFile) -> Result<NeoData> {
        let info = source.variant_info.clone();
        let mut chromosomes = [String::new(), String::new()];
        let mut positions = [0, 0];
        let mut orientations: [i8; 2] = [0, 0];
        if info.contains(';') {
            let streams: Vec<&str> = hmf_common::java_text::split(&info, ';');
            for fs in [FS_UP, FS_DOWN] {
                let parts: Vec<&str> = hmf_common::java_text::split(streams[fs], ':');
                chromosomes[fs] = parts[0].to_string();
                positions[fs] = parts[1].parse().with_context(|| format!("invalid variant info {info}"))?;
                orientations[fs] = parts[2].parse().with_context(|| format!("invalid variant info {info}"))?;
            }
        } else {
            let parts: Vec<&str> = hmf_common::java_text::split(&info, ':');
            chromosomes = [parts[0].to_string(), parts[0].to_string()];
            let p: i32 = parts[1].parse().with_context(|| format!("invalid variant info {info}"))?;
            positions = [p, p];
        }
        let transcripts = [
            hmf_common::java_text::split(&source.transcripts[0], ';').into_iter().map(|s| s.to_string()).collect(),
            hmf_common::java_text::split(&source.transcripts[1], ';').into_iter().map(|s| s.to_string()).collect(),
        ];
        let mut coords: [Vec<[i32; 2]>; 2] = [Vec::new(), Vec::new()];
        for fs in [FS_UP, FS_DOWN] {
            if source.coding_base_cigars[fs].is_empty() {
                continue;
            }
            let c = parse_cigar(&source.coding_base_cigars[fs]).unwrap_or_default();
            coords[fs] = generate_mapped_coords(&c, source.coding_base_positions[fs][0]);
        }
        Ok(NeoData { source, chromosomes, positions, orientations, transcripts, coding_base_coords: coords, support: Support::default() })
    }
    fn is_fusion(&self) -> bool {
        self.source.variant_type.is_fusion()
    }
    fn is_point_mutation(&self) -> bool {
        self.source.variant_type.is_point_mutation()
    }
    fn single_gene(&self) -> bool {
        self.source.gene_ids[0] == self.source.gene_ids[1]
    }
    fn single_chromosome(&self) -> bool {
        self.chromosomes[0] == self.chromosomes[1]
    }
    fn coding_base_range(&self, fs: usize) -> [i32; 2] {
        self.source.coding_base_positions[fs]
    }
    fn full_coding_bases(&self, perspective: usize) -> String {
        let cb = &self.source.coding_bases;
        let rc = |s: &str| String::from_utf8(hmf_common::codon::reverse_complement_bases(s.as_bytes())).unwrap_or_default();
        if !self.is_fusion() || (self.orientations[0] != self.orientations[1] && self.single_chromosome()) {
            return if self.orientations[0] == ORIENT_FWD { format!("{}{}", cb[0], cb[1]) } else { format!("{}{}", cb[1], cb[0]) };
        }
        if perspective == FS_UP {
            if self.orientations[0] == ORIENT_FWD { format!("{}{}", cb[0], rc(&cb[1])) } else { format!("{}{}", rc(&cb[1]), cb[0]) }
        } else if self.orientations[1] == ORIENT_FWD {
            format!("{}{}", cb[1], rc(&cb[0]))
        } else {
            format!("{}{}", rc(&cb[0]), cb[1])
        }
    }
    fn is_deletion_fusion(&self) -> bool {
        if !self.is_fusion() || !self.single_chromosome() {
            return false;
        }
        if self.positions[0] < self.positions[1] {
            self.orientations[0] == ORIENT_FWD && self.orientations[1] == ORIENT_REV
        } else {
            self.orientations[0] == ORIENT_REV && self.orientations[1] == ORIENT_FWD
        }
    }
    fn fusion_soft_clipped_bases(&self, perspective: usize) -> String {
        let other = &self.source.coding_bases[switch_index(perspective)];
        if self.orientations[0] != self.orientations[1] {
            other.clone()
        } else {
            String::from_utf8(hmf_common::codon::reverse_complement_bases(other.as_bytes())).unwrap_or_default()
        }
    }
}

fn calc_base_overlap(r: &BaseRegion, c: &[i32; 2]) -> i32 {
    let s = r.start.max(c[0]);
    let e = r.end.min(c[1]);
    if s <= e { e - s + 1 } else { 0 }
}

fn calc_coordinates_overlap(regions: &[BaseRegion], coords: &[[i32; 2]]) -> i32 {
    let mut n = 0;
    for r in regions {
        for c in coords {
            n += calc_base_overlap(r, c);
        }
    }
    n
}

fn calc_base_match(b1: &str, b2: &str) -> usize {
    if b1 == b2 || b1.contains(b2) || b2.contains(b1) {
        return EXACT_MATCH;
    }
    let mut mismatches = 0;
    let n = b1.len().min(b2.len());
    let (x, y) = (b1.as_bytes(), b2.as_bytes());
    for i in 0..n {
        if x[i] != y[i] {
            mismatches += 1;
            if mismatches > MAX_BASE_MISMATCH {
                return MISMATCH;
            }
        }
    }
    PARTIAL_MATCH
}

fn compare_coding_bases(read: &Read, neo_bases: &str, neo_coords: &[[i32; 2]], pos_start: i32, pos_end: i32) -> usize {
    let mut idx: i32 = 0;
    let mut read_bases = String::new();
    if read.is_soft_clipped_no_region_match(SE_START) {
        idx += read.left_clip_length();
    }
    let rb = read.read_bases();
    for c in &read.mapped_coords.alignments {
        if pos_start > c.end {
            idx += c.end - c.start + 1;
            continue;
        }
        let mut pos = if position_within(pos_start, c.start, c.end) {
            idx += pos_start - c.start;
            pos_start
        } else {
            c.start
        };
        while pos <= c.end {
            if pos > pos_end || idx >= read.base_length() {
                break;
            }
            read_bases.push(rb[idx as usize] as char);
            idx += 1;
            pos += 1;
        }
        if pos > pos_end {
            break;
        }
    }
    let mut nidx: i32 = 0;
    let mut neo = String::new();
    let nb = neo_bases.as_bytes();
    for c in neo_coords {
        if pos_start > c[1] {
            nidx += c[1] - c[0] + 1;
            continue;
        }
        let mut pos = if position_within(pos_start, c[0], c[1]) {
            nidx += pos_start - c[0];
            pos_start
        } else {
            c[0]
        };
        while pos <= c[1] {
            if pos > pos_end || nidx >= nb.len() as i32 {
                break;
            }
            neo.push(nb[nidx as usize] as char);
            nidx += 1;
            pos += 1;
        }
        if pos > pos_end {
            break;
        }
    }
    calc_base_match(&read_bases, &neo)
}

/// NeoFragmentMatcher.findFusionSupport
fn find_fusion_support(neo: &NeoData, stream: usize, read: &Read) -> Support {
    let mut s = Support::default();
    let junction_side = if neo.orientations[stream] == ORIENT_FWD { SE_END } else { SE_START };
    let range = neo.coding_base_range(stream);
    let within = (junction_side == SE_START && read.pos_start >= range[0]) || (junction_side == SE_END && read.pos_end <= range[1]);
    let coords = &read.mapped_coords.alignments;
    if within {
        let bases = &neo.source.coding_bases[stream];
        if calc_coordinates_overlap(coords, &neo.coding_base_coords[stream]) < MIN_BASE_OVERLAP {
            return s;
        }
        let max_start = read.pos_start.max(range[0]);
        let min_end = read.pos_end.min(range[1]);
        let mut level = compare_coding_bases(read, bases, &neo.coding_base_coords[stream], max_start, min_end);
        if level == MISMATCH {
            return s;
        }
        if read.is_soft_clipped_no_region_match(junction_side) {
            let post = neo.fusion_soft_clipped_bases(stream);
            let sc = if junction_side == SE_START { read.left_clip_length() } else { read.right_clip_length() };
            if sc >= MIN_BASE_OVERLAP / 2 {
                let rb = read.read_bases();
                let clipped = if junction_side == SE_START { &rb[..sc as usize] } else { &rb[rb.len() - sc as usize..] };
                level = calc_base_match(&post, std::str::from_utf8(clipped).unwrap_or(""));
                if level == MISMATCH {
                    return s;
                }
                s.novel[level] += 1;
            }
        } else if stream == FS_UP {
            s.up[level] += 1;
        } else {
            s.down[level] += 1;
        }
    } else {
        if !neo.is_deletion_fusion() {
            return s;
        }
        let bases = neo.full_coding_bases(stream);
        let junction = if neo.orientations[FS_UP] == ORIENT_FWD {
            [neo.source.coding_base_positions[FS_UP][1], neo.source.coding_base_positions[FS_DOWN][0]]
        } else {
            [neo.source.coding_base_positions[FS_DOWN][1], neo.source.coding_base_positions[FS_UP][0]]
        };
        let mut supports_split = false;
        for i in 0..coords.len().saturating_sub(1) {
            if coords[i].end == junction[0] && coords[i + 1].start == junction[1] {
                supports_split = true;
                break;
            }
        }
        if !supports_split {
            return s;
        }
        if calc_coordinates_overlap(coords, &neo.coding_base_coords[stream]) < MIN_BASE_OVERLAP {
            return s;
        }
        let max_start = read.pos_start.max(range[0]);
        let min_end = read.pos_end.min(range[1]);
        let level = compare_coding_bases(read, &bases, &neo.coding_base_coords[stream], max_start, min_end);
        if level == MISMATCH {
            return s;
        }
        if position_within(neo.positions[FS_UP], read.pos_start, read.pos_end) {
            s.novel[level] += 1;
        } else if stream == FS_UP {
            s.up[level] += 1;
        } else {
            s.down[level] += 1;
        }
    }
    s
}

fn check_base_coverage(neo: &mut NeoData, reads: &[&Read]) {
    for fs in [FS_UP, FS_DOWN] {
        if neo.is_point_mutation() && fs == FS_DOWN {
            neo.support.ref_base_depth[FS_DOWN] = neo.support.ref_base_depth[FS_UP];
            break;
        }
        let ref_base = if neo.is_point_mutation() {
            neo.positions[fs]
        } else {
            neo.coding_base_range(fs)[if neo.orientations[fs] == ORIENT_FWD { SE_END } else { SE_START }]
        };
        let covers = reads.iter().any(|r| r.chromosome == neo.chromosomes[fs] && r.mapped_coords.alignments.iter().any(|x| position_within(ref_base, x.start, x.end)));
        if covers {
            neo.support.ref_base_depth[fs] += 1;
        }
    }
}

/// a read and the gene collection (0 up, 1 down) whose regions it was matched against
struct NeoRead {
    read: Read,
    gc: usize,
}

struct ReadGroup {
    reads: Vec<usize>,
    complete: bool,
}

fn group_complete(reads: &[usize], store: &[NeoRead]) -> bool {
    let mut supp = 0;
    let mut non_supp = 0;
    let mut exp_supp = 0;
    let mut exp_non_supp = 1;
    for i in reads {
        let r = &store[*i].read;
        if r.is_read_paired() && !r.is_mate_unmapped() {
            exp_non_supp = 2;
        }
        if r.is_supplementary() {
            supp += 1;
        } else {
            non_supp += 1;
            if r.has_supp_alignment() {
                exp_supp += 1;
            }
        }
    }
    exp_non_supp == non_supp && exp_supp == supp
}

fn process_fragment_reads(neo: &mut NeoData, group: &[usize], store: &[NeoRead], gcs: &[Option<GeneCollection>; 2]) {
    let reads: Vec<&Read> = group.iter().map(|i| &store[*i].read).collect();
    if reads.iter().any(|r| r.is_duplicate()) {
        return;
    }
    check_base_coverage(neo, &reads);
    if neo.is_fusion() && !reads.iter().any(|r| r.has_supp_alignment()) && !reads.iter().any(|r| r.has_split) {
        return;
    }
    let mut supported = [false, false];
    for fs in [FS_UP, FS_DOWN] {
        if neo.single_gene() && fs == FS_DOWN {
            supported[FS_DOWN] = supported[FS_UP];
            break;
        }
        for name in &neo.transcripts[fs] {
            let mut has_trans = false;
            let mut has_support = true;
            for i in group {
                let nr = &store[*i];
                let Some(gc) = &gcs[nr.gc] else { continue };
                for (ri, mt) in nr.read.mapped_regions.entries() {
                    let region = &gc.exon_regions[ri];
                    if !region.trans_exon_refs.iter().any(|x| x.trans_name == name.as_str()) {
                        continue;
                    }
                    has_trans = true;
                    if !mt.valid_exon_match() {
                        has_support = false;
                    }
                }
            }
            if has_trans && has_support {
                supported[fs] = true;
                break;
            }
        }
    }
    if !supported[0] || !supported[1] {
        return;
    }
    let mut group_support = Support::default();
    for fs in [FS_UP, FS_DOWN] {
        let range = neo.coding_base_range(fs);
        for r in &reads {
            if r.chromosome != neo.chromosomes[fs] {
                continue;
            }
            if !r.mapped_coords.without_inferred().iter().any(|x| positions_overlap(range[0], range[1], x.start, x.end)) {
                continue;
            }
            let s = find_fusion_support(neo, fs, r);
            group_support.set_max(&s);
        }
    }
    neo.support.combine(&group_support);
}

/// NeoEpitopeReader.calcFragmentSupport
pub fn calc_fragment_support(config: &IsofoxConfig, cache: &'static EnsemblDataCache) -> Result<()> {
    let neo_dir = config.neo_dir.clone().unwrap_or_default();
    let tumor = config.dna_tumor_id.clone().unwrap_or_else(|| "null".to_string());
    let filename = format!("{}{}.neo.neo_data.tsv", hmf_common::config_builder::check_add_dir_separator(&neo_dir), tumor);
    let mut neos: Vec<NeoData> = Vec::new();
    match read_neo_epitope_file(&filename) {
        Ok(list) => {
            for n in list {
                neos.push(NeoData::new(n)?);
            }
            hmf_common::info!("sample({}) loaded {} neo-epitopes from file: {}", config.sample_id, neos.len(), filename);
        }
        Err(e) => hmf_common::error!("failed to read neo-epitope file({}): {}", filename, e),
    }
    let out_file = format!("{}{}.isf.neoepitope.tsv", config.output_dir, config.sample_id);
    let mut writer = crate::run::Writer::create(&out_file)?;
    writer.line(&rna_neo_epitope_header())?;

    let mut reader = match &config.bam_file {
        Some(b) => Some(BamReader::open(b, config.ref_genome_file.as_deref())?),
        None => None,
    };
    let mut refg = RefSource::open(config.ref_genome_file.as_deref().unwrap_or(""))?;
    let filter = SliceFilter { min_map_quality: config.single_map_quality, keep_duplicates: true, keep_supplementaries: true, keep_secondaries: false };
    let restricted = &config.filters.restricted_gene_ids;

    for neo in &mut neos {
        if !restricted.is_empty() && !(restricted.contains(&neo.source.gene_ids[0]) && restricted.contains(&neo.source.gene_ids[1])) {
            continue;
        }
        if neo.is_point_mutation() {
            continue;
        }
        let mut store: Vec<NeoRead> = Vec::new();
        let mut groups: SeqHashMap<String, ReadGroup> = SeqHashMap::new();
        let mut gcs: [Option<GeneCollection>; 2] = [None, None];
        for fs in [FS_UP, FS_DOWN] {
            let gene_id = neo.source.gene_ids[fs].clone();
            let Some(gene) = cache.gene_data_by_id(&gene_id) else {
                // Java logs the missing gene and then fails on the next read
                hmf_common::error!("gene({}) not found", gene_id);
                anyhow::bail!("gene({gene_id}) not found");
            };
            let list = create_gene_read_data(&[gene], cache);
            if list.is_empty() {
                anyhow::bail!("gene({gene_id}) has no transcripts");
            }
            let mut gc = GeneCollection::new(0, list);
            if config.ref_genome_file.is_some() {
                for r in &mut gc.exon_regions {
                    let b = refg.bases(gc.chromosome, r.start, r.end);
                    r.set_ref_bases(b);
                }
            }
            let range = neo.source.coding_base_positions[fs];
            let chromosome = neo.chromosomes[fs].clone();
            let mut processed_complete: Vec<Vec<usize>> = Vec::new();
            if let Some(rd) = reader.as_mut() {
                rd.slice(&chromosome, range[0], range[1], &filter, |rec: &bam::Record| {
                    let mut read = to_read(rec);
                    let overlapping = find_overlapping_regions(&gc.exon_regions, None, &read);
                    process_overlapping_regions(&mut read, &gc.exon_regions, &overlapping);
                    gc.set_read_gene_collections(&mut read, gc.region_bounds);
                    let id = read.id.to_string();
                    store.push(NeoRead { read, gc: fs });
                    let ri = store.len() - 1;
                    match groups.get_mut(&id) {
                        None => {
                            let complete = group_complete(&[ri], &store);
                            groups.insert(id, ReadGroup { reads: vec![ri], complete });
                        }
                        Some(g) => {
                            g.reads.push(ri);
                            g.complete = group_complete(&g.reads, &store);
                            if g.complete {
                                processed_complete.push(g.reads.clone());
                                groups.remove(&id);
                            }
                        }
                    }
                    true
                })?;
            }
            gcs[fs] = Some(gc);
            for g in processed_complete {
                process_fragment_reads(neo, &g, &store, &gcs);
            }
        }
        let remaining: Vec<Vec<usize>> = groups.values().map(|g| g.reads.clone()).collect();
        for g in remaining {
            process_fragment_reads(neo, &g, &store, &gcs);
        }
        writer.line(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            neo.source.id,
            neo.source.variant_type.name(),
            neo.source.variant_info,
            neo.support.novel[EXACT_MATCH],
            neo.support.ref_base_depth[FS_UP],
            neo.support.ref_base_depth[FS_DOWN]
        ))?;
    }
    writer.close()
}
