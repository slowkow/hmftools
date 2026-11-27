//! CiderReadScreener: which reads are V or J anchor reads (VJReadCandidate), by their alignment to an
//! anchor location, by a BLOSUM anchor search in the soft clip of reads at a constant or D region, or in
//! unmapped reads whose mate maps near an anchor; then the mates of those reads. Also CiderUtils'
//! adapter trimming.

use anyhow::{bail, Result};
use hmf_common::codon::{amino_acid_from_bases, reverse_complement};

use crate::blosum::{similarity_score, AnchorSearcher, SearchMode};
use crate::compat;
use crate::genes::{GeneDatastore, GenomicLocation, IgTcrConstantDiversityRegion, Strand, Vj, VjAnchorGenomeLocation, VjGeneType};
use crate::read::{left_soft_clip, right_soft_clip, CiderRead, ReadRef};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MatchMethod {
    Align,
    Blosum,
}

/// VJReadCandidate; the templates are indices into the datastore's templates
#[derive(Clone, Debug)]
pub struct VjReadCandidate {
    pub read: ReadRef,
    pub templates: Vec<usize>,
    pub gene_type: VjGeneType,
    pub template_anchor_sequence: String,
    pub match_method: MatchMethod,
    pub use_reverse_complement: bool,
    /// after the reverse complement if used; may be negative
    pub anchor_offset_start: i32,
    /// after the reverse complement if used; may be past the end
    pub anchor_offset_end: i32,
    pub left_soft_clip: i32,
    pub right_soft_clip: i32,
    pub similarity_score: i32,
}

impl VjReadCandidate {
    pub fn read_length(&self) -> i32 {
        self.read.read_length()
    }
    /// the read in the order the gene is transcribed
    pub fn read_sequence(&self) -> Vec<u8> {
        let s = self.read.read_string();
        if self.use_reverse_complement { reverse_complement(s) } else { s.to_vec() }
    }
    pub fn anchor_sequence(&self) -> Vec<u8> {
        let len = self.read_length();
        if self.anchor_offset_start >= len || self.anchor_offset_end <= 0 {
            return Vec::new();
        }
        let seq = self.read_sequence();
        seq[self.anchor_offset_start.max(0) as usize..self.anchor_offset_end.min(len) as usize].to_vec()
    }
#[allow(dead_code)]
    pub fn anchor_aa(&self) -> String {
        amino_acid_from_bases(&self.anchor_sequence())
    }
}

/// SamRecordUtils.getFivePrimeUnclippedPosition(readStart, cigarStr, forwardStrand) and
/// CigarUtils.getReadBoundaryPosition(readStart, cigarStr, getReadStart, includeSoftClipped = true)
pub fn read_boundary_position(read_start: i32, cigar: &str, get_read_start: bool) -> i32 {
    let mut current = read_start;
    let mut element_length: i32 = 0;
    for c in cigar.bytes() {
        let add = matches!(c, b'D' | b'M' | b'S' | b'N');
        if add {
            if get_read_start {
                return if c == b'S' { read_start - element_length } else { read_start };
            }
            if c == b'S' && read_start == current {
                // a left clip, ignored for the end position
            } else {
                current = current.wrapping_add(element_length);
            }
            element_length = 0;
            continue;
        }
        if c.is_ascii_digit() {
            element_length = element_length.wrapping_mul(10).wrapping_add((c - b'0') as i32);
        } else {
            element_length = 0;
        }
    }
    current - 1
}

/// CiderUtils.calculateAdapterSequenceTrim: the bases to trim from the left and right of a read that runs
/// into adapter (a fragment shorter than the read).
///
/// CI2: Java takes the read's 5' position with its negative strand flag as the forward strand flag, so it
/// gets the read's 3' end, and the insert length is that of the mate's 5' end to the read's 3' end.
/// Corrected: the read's 5' end.
pub fn adapter_sequence_trim(read: &CiderRead) -> (i32, i32) {
    if !read.unmapped() && !read.mate_unmapped() && read.tid >= 0 && read.contig == read.mate_contig && read.negative_strand() != read.mate_negative_strand() {
        if let Some(mate_cigar) = &read.mate_cigar {
            let cigar = read.cigar_string();
            let forward_flag = if compat::java("CI2") { read.negative_strand() } else { !read.negative_strand() };
            let five_prime = read_boundary_position(read.alignment_start, &cigar, forward_flag);
            let three_prime =
                if read.negative_strand() { read.alignment_start - left_soft_clip(&read.cigar) } else { read.alignment_end() + right_soft_clip(&read.cigar) };
            let mate_five_prime = read_boundary_position(read.mate_alignment_start, mate_cigar, !read.mate_negative_strand());
            let insert_length = (mate_five_prime - five_prime).abs();
            let adapter_bases = if read.negative_strand() { (mate_five_prime - three_prime).max(0) } else { (three_prime - mate_five_prime).max(0) };
            if insert_length < read.read_length() && adapter_bases > 0 {
                return if read.negative_strand() { (adapter_bases, 0) } else { (0, adapter_bases) };
            }
        }
    }
    (0, 0)
}

/// a mapped region: chromosome, start, end (GenomeRegions.create)
pub struct Mapped<'a> {
    pub chromosome: &'a str,
    pub start: i32,
    pub end: i32,
}

pub fn mapped_region(read: &CiderRead) -> Mapped<'_> {
    Mapped { chromosome: &read.contig, start: read.alignment_start, end: read.alignment_end() }
}

/// CiderReadScreener.isMappedToAnchorLocation
pub fn is_mapped_to_anchor_location(mapped: &Mapped, loc: &VjAnchorGenomeLocation, tolerance: i32) -> bool {
    if loc.location.chromosome != mapped.chromosome {
        return false;
    }
    loc.start() <= mapped.end + tolerance && mapped.start - tolerance <= loc.end()
}

/// CiderReadScreener.extrapolateReadOffsetAtRefPosition: the 0-based read offset at a reference position,
/// extrapolated from the closest alignment block (may be negative or past the read).
///
/// CI4: Java measures the distance to a block that ends before the position from the block's start, not
/// its end, so a position in a deletion or splice gap goes to the next block though it is closer to the
/// block before. Corrected: from the block's end.
pub fn extrapolate_read_offset_at_ref_position(read: &CiderRead, ref_pos: i32) -> i32 {
    let mut closest = i32::MAX;
    let mut read_offset = 0;
    for b in read.alignment_blocks() {
        let block_end = b.reference_start + b.length - 1;
        let block_start = b.reference_start;
        let distance = if ref_pos < block_start {
            block_start - ref_pos
        } else if ref_pos > block_end {
            if compat::java("CI4") { ref_pos - block_start } else { ref_pos - block_end }
        } else {
            0
        };
        if distance < closest {
            closest = distance;
            read_offset = b.read_start + ref_pos - block_start - 1;
        }
    }
    read_offset
}

/// CiderReadScreener.extrapolateAnchorReadRange: (start, end exclusive) of the anchor in the read
pub fn extrapolate_anchor_read_range(read: &CiderRead, loc: &VjAnchorGenomeLocation) -> Result<(i32, i32)> {
    let ref_pos = loc.anchor_boundary_reference_position();
    let offset = extrapolate_read_offset_at_ref_position(read, ref_pos);
    if ref_pos == loc.start() {
        Ok((offset, offset + loc.base_length()))
    } else if ref_pos == loc.end() {
        Ok((offset - loc.base_length() + 1, offset + 1))
    } else {
        bail!("anchor end reference pos must be either anchor start or anchor end")
    }
}

/// isUnamppedReadRelevantToAnchorLoc: an unmapped read whose mate maps upstream of a V anchor or downstream
/// of a J anchor, pointing at it, within the fragment length
pub fn is_unmapped_read_relevant_to_anchor_loc(read: &CiderRead, loc: &VjAnchorGenomeLocation, max_fragment_length: i32) -> bool {
    if !read.paired() || read.mate_unmapped_flag() {
        return false;
    }
    if loc.location.chromosome.as_str() != &*read.mate_contig {
        return false;
    }
    let mate_forward = !read.mate_negative_strand();
    let mate_start = read.mate_alignment_start;
    let inferred_mate_end = mate_start + read.read_length();
    if (loc.vj() == Vj::V && loc.strand() == Strand::Forward) || (loc.vj() == Vj::J && loc.strand() == Strand::Reverse) {
        if !mate_forward || mate_start > loc.end() {
            return false;
        }
        (loc.start() - inferred_mate_end) < max_fragment_length
    } else {
        if mate_forward || inferred_mate_end < loc.start() {
            return false;
        }
        (mate_start - loc.end()) < max_fragment_length
    }
}

/// isUnamppedReadRelevantToConstantRegion
pub fn is_unmapped_read_relevant_to_constant_region(read: &CiderRead, loc: &GenomicLocation, max_fragment_length: i32) -> bool {
    if !read.paired() || read.mate_unmapped_flag() {
        return false;
    }
    if loc.chromosome.as_str() != &*read.mate_contig {
        return false;
    }
    let mate_start = read.mate_alignment_start;
    let inferred_mate_end = mate_start + read.read_length();
    if (mate_start - loc.pos_end) >= max_fragment_length || (loc.pos_start - inferred_mate_end) >= max_fragment_length {
        return false;
    }
    let mate_strand = if read.mate_negative_strand() { Strand::Reverse } else { Strand::Forward };
    loc.strand != mate_strand
}

/// The outcome of screening one record
pub enum Screened {
    Candidate(VjReadCandidate),
    /// mapped, not a candidate: its mate may make it one
    UnmatchedMapped,
    Unmatched,
}

/// The anchor locations and constant regions of each chromosome, for the per-read loops: Java tries them
/// all in list order and takes the first match; only those on the read's chromosome can match, and among
/// those the list order is kept.
struct ChromosomeIndex {
    /// per chromosome: (start, end, index in the anchor location list), by start
    anchors: std::collections::HashMap<String, Vec<(i32, i32, usize)>>,
    /// the longest anchor location
    max_anchor_length: i32,
    /// per chromosome: the anchor location indices in list order
    anchors_in_order: std::collections::HashMap<String, Vec<usize>>,
    /// per chromosome: the constant region indices in list order
    regions_in_order: std::collections::HashMap<String, Vec<usize>>,
}

impl ChromosomeIndex {
    fn new(ds: &GeneDatastore) -> Self {
        let mut anchors: std::collections::HashMap<String, Vec<(i32, i32, usize)>> = std::collections::HashMap::new();
        let mut anchors_in_order: std::collections::HashMap<String, Vec<usize>> = std::collections::HashMap::new();
        let mut max_anchor_length = 0;
        for (i, l) in ds.anchor_locations.iter().enumerate() {
            anchors.entry(l.location.chromosome.clone()).or_default().push((l.start(), l.end(), i));
            anchors_in_order.entry(l.location.chromosome.clone()).or_default().push(i);
            max_anchor_length = max_anchor_length.max(l.end() - l.start() + 1);
        }
        for v in anchors.values_mut() {
            v.sort();
        }
        let mut regions_in_order: std::collections::HashMap<String, Vec<usize>> = std::collections::HashMap::new();
        for (i, r) in ds.constant_diversity_regions.iter().enumerate() {
            regions_in_order.entry(r.location.chromosome.clone()).or_default().push(i);
        }
        ChromosomeIndex { anchors, max_anchor_length, anchors_in_order, regions_in_order }
    }

    /// the anchor locations overlapping [start, end] on the chromosome, in list order
    fn overlapping_anchors(&self, chromosome: &str, start: i32, end: i32, out: &mut Vec<usize>) {
        out.clear();
        let Some(v) = self.anchors.get(chromosome) else { return };
        let from = v.partition_point(|a| a.0 < start - self.max_anchor_length);
        for a in &v[from..] {
            if a.0 > end {
                break;
            }
            if a.0 <= end && start <= a.1 {
                out.push(a.2);
            }
        }
        out.sort_unstable();
    }
}

pub struct CiderReadScreener<'a> {
    ds: &'a GeneDatastore,
    searcher: &'a dyn AnchorSearcher,
    min_read_cdr3_overlap: i32,
    max_fragment_length: i32,
    index: ChromosomeIndex,
}

impl<'a> CiderReadScreener<'a> {
    pub fn new(ds: &'a GeneDatastore, searcher: &'a dyn AnchorSearcher, min_read_cdr3_overlap: i32, max_fragment_length: i32) -> Self {
        CiderReadScreener { ds, searcher, min_read_cdr3_overlap, max_fragment_length, index: ChromosomeIndex::new(ds) }
    }

    /// tryMatchRead; `read` is the shared record, for the candidate
    pub fn screen(&self, read: &ReadRef) -> Result<Screened> {
        let fragment_mapped = read.tid != -1;
        let read_mapped = !read.unmapped();
        let candidate = if fragment_mapped {
            if read_mapped { self.try_match_by_alignment(read)? } else { self.try_match_unmapped_read_by_blosum(read)? }
        } else {
            self.try_match_by_blosum(read)?
        };
        Ok(match candidate {
            Some(c) => Screened::Candidate(c),
            None if read_mapped => Screened::UnmatchedMapped,
            None => Screened::Unmatched,
        })
    }

    fn try_match_by_alignment(&self, read: &ReadRef) -> Result<Option<VjReadCandidate>> {
        let mapped = mapped_region(read);
        // tryMatchToAnchorLocationExactly: the anchor locations the read overlaps, in list order
        let mut overlapping = Vec::new();
        self.index.overlapping_anchors(mapped.chromosome, mapped.start, mapped.end, &mut overlapping);
        for &i in &overlapping {
            if let Some(c) = self.try_match_to_anchor_location(read, &mapped, &self.ds.anchor_locations[i], 0)? {
                return Ok(Some(c));
            }
        }
        if left_soft_clip(&read.cigar) != 0 || right_soft_clip(&read.cigar) != 0 {
            if let Some(regions) = self.index.regions_in_order.get(mapped.chromosome) {
                for &i in regions {
                    if let Some(c) = self.try_match_from_constant_region(read, &mapped, &self.ds.constant_diversity_regions[i])? {
                        return Ok(Some(c));
                    }
                }
            }
        }
        Ok(None)
    }

    /// tryMatchToAnchorLocation
    pub fn try_match_to_anchor_location(
        &self,
        read: &ReadRef,
        mapped: &Mapped,
        loc: &VjAnchorGenomeLocation,
        tolerance: i32,
    ) -> Result<Option<VjReadCandidate>> {
        if !is_mapped_to_anchor_location(mapped, loc, tolerance) {
            return Ok(None);
        }
        if loc.base_length() != 30 {
            bail!("unexpected anchor length");
        }
        let (mut start, mut end) = extrapolate_anchor_read_range(read, loc)?;
        let trim = read.adapter_trim;
        let near = (loc.anchor_boundary_side() == 1 && read.read_length() - trim.1 - end >= self.min_read_cdr3_overlap - tolerance)
            || (loc.anchor_boundary_side() == -1 && start - trim.0 >= self.min_read_cdr3_overlap - tolerance);
        if near {
            if loc.strand() == Strand::Reverse {
                let rev_start = read.read_length() - end;
                end = read.read_length() - start;
                start = rev_start;
            }
            let genes = self.ds.by_location(&loc.location);
            if !genes.is_empty() {
                return self.create_candidate(read, genes, MatchMethod::Align, loc.strand() == Strand::Reverse, start, end).map(Some);
            }
        }
        Ok(None)
    }

    /// tryMatchUnmappedReadByBlosum: an unmapped read whose mate maps near an anchor or a constant region;
    /// the anchor search is for that locus only
    fn try_match_unmapped_read_by_blosum(&self, read: &ReadRef) -> Result<Option<VjReadCandidate>> {
        if !read.paired() || !read.unmapped() || read.mate_unmapped_flag() {
            bail!("java.lang.IllegalArgumentException: Failed requirement (unmapped read with a mapped mate)");
        }
        // the first anchor location in list order on the mate's chromosome that is relevant
        let anchor = self
            .index
            .anchors_in_order
            .get(&*read.mate_contig)
            .and_then(|v| v.iter().map(|&i| &self.ds.anchor_locations[i]).find(|l| is_unmapped_read_relevant_to_anchor_loc(read, l, self.max_fragment_length)));
        let mut types: Vec<VjGeneType> = Vec::new();
        match anchor {
            Some(l) => {
                types.push(l.gene_type);
                types.extend_from_slice(l.gene_type.paired());
            }
            None => {
                let region = self.index.regions_in_order.get(&*read.mate_contig).and_then(|v| {
                    v.iter().map(|&i| &self.ds.constant_diversity_regions[i]).find(|r| is_unmapped_read_relevant_to_constant_region(read, &r.location, self.max_fragment_length))
                });
                match region {
                    Some(r) => types.extend_from_slice(r.corresponding_vj()),
                    None => return Ok(None),
                }
            }
        }
        for reverse in [true, false] {
            let s = if reverse { reverse_complement(read.read_string()) } else { read.read_string().to_vec() };
            let m = self.searcher.search(&s, &types, SearchMode::DisallowNegSimilarity, 0, read.read_length())?;
            if let Some(m) = m {
                if m.similarity_score > 0 && !m.template_genes.is_empty() {
                    return self.create_candidate(read, &m.template_genes, MatchMethod::Blosum, reverse, m.anchor_start, m.anchor_end).map(Some);
                }
            }
        }
        Ok(None)
    }

    /// tryMatchByBlosum: a read of an unmapped pair, both strands, all gene types
    fn try_match_by_blosum(&self, read: &ReadRef) -> Result<Option<VjReadCandidate>> {
        for reverse in [false, true] {
            let s = if reverse { reverse_complement(read.read_string()) } else { read.read_string().to_vec() };
            if let Some(m) = self.searcher.search_all(&s, SearchMode::DisallowNegSimilarity)? {
                if m.similarity_score > 0 && !m.template_genes.is_empty() {
                    return self.create_candidate(read, &m.template_genes, MatchMethod::Blosum, reverse, m.anchor_start, m.anchor_end).map(Some);
                }
            }
        }
        Ok(None)
    }

    /// tryMatchFromConstantRegion: a read at a constant or D region with a soft clip, which may hold a J
    /// (or V) anchor
    pub fn try_match_from_constant_region(&self, read: &ReadRef, mapped: &Mapped, region: &IgTcrConstantDiversityRegion) -> Result<Option<VjReadCandidate>> {
        if !region.location.in_primary_assembly() {
            return Ok(None);
        }
        let read_length = read.read_length();
        let loc = &region.location;
        if loc.pos_start - read_length < mapped.end && loc.pos_end + read_length > mapped.start && loc.chromosome == mapped.chromosome {
            let m = if loc.strand == Strand::Forward {
                let clip = left_soft_clip(&read.cigar);
                if clip == 0 {
                    return Ok(None);
                }
                self.searcher.search(read.read_string(), region.corresponding_vj(), SearchMode::DisallowNegSimilarity, 0, clip)?
            } else {
                let clip = right_soft_clip(&read.cigar);
                if clip == 0 {
                    return Ok(None);
                }
                let rc = reverse_complement(read.read_string());
                self.searcher.search(&rc, region.corresponding_vj(), SearchMode::DisallowNegSimilarity, 0, clip)?
            };
            if let Some(m) = m {
                if m.similarity_score > 0 && !m.template_genes.is_empty() {
                    return self
                        .create_candidate(read, &m.template_genes, MatchMethod::Blosum, loc.strand == Strand::Reverse, m.anchor_start, m.anchor_end)
                        .map(Some);
                }
            }
        }
        Ok(None)
    }

    /// createVjReadCandidate
    fn create_candidate(&self, read: &ReadRef, templates: &[usize], method: MatchMethod, rev: bool, start: i32, end: i32) -> Result<VjReadCandidate> {
        let first = self.ds.template(templates[0]);
        let gene_type = first.gene_type;
        if templates.iter().any(|&t| self.ds.template(t).gene_type != gene_type) {
            crate::error!("multiple gene types found in same match");
            bail!("multiple gene types found in same match");
        }
        let mut c = VjReadCandidate {
            read: read.clone(),
            templates: templates.to_vec(),
            gene_type,
            template_anchor_sequence: first.anchor_sequence.clone(),
            match_method: method,
            use_reverse_complement: rev,
            anchor_offset_start: start,
            anchor_offset_end: end,
            left_soft_clip: left_soft_clip(&read.cigar),
            right_soft_clip: right_soft_clip(&read.cigar),
            similarity_score: i32::MIN,
        };
        if method == MatchMethod::Blosum {
            c.similarity_score = similarity_score(gene_type.vj(), c.template_anchor_sequence.as_bytes(), &c.anchor_sequence())?;
            if c.similarity_score <= 0 {
                crate::error!("blosum match with -ve similarity score");
                bail!("blosum match with -ve similarity score");
            }
        }
        Ok(c)
    }

    /// processCandidateMates: each unmatched mapped read whose fragment has candidates is matched, with
    /// the fragment length as tolerance, to the anchor locations of those candidates; the closest wins
    pub fn process_candidate_mates(&self, candidates: &[VjReadCandidate], unmatched_mapped: &[ReadRef]) -> Result<Vec<VjReadCandidate>> {
        let mut by_fragment: std::collections::HashMap<&str, Vec<VjAnchorGenomeLocation>> = std::collections::HashMap::new();
        for c in candidates {
            let entry = by_fragment.entry(c.read.name.as_str()).or_default();
            for &t in &c.templates {
                let template = self.ds.template(t);
                if let Some(loc) = &template.anchor_location {
                    entry.push(VjAnchorGenomeLocation { gene_type: template.gene_type, location: loc.clone() });
                }
            }
        }
        let mut mates = Vec::new();
        for read in unmatched_mapped {
            let Some(locations) = by_fragment.get(read.name.as_str()) else { continue };
            let mapped = mapped_region(read);
            let mut best: Option<(VjReadCandidate, i32)> = None;
            for loc in locations {
                if let Some(c) = self.try_match_to_anchor_location(read, &mapped, loc, self.max_fragment_length)? {
                    let distance = c.anchor_offset_start.abs().min(c.anchor_offset_end.abs());
                    if best.as_ref().is_none_or(|b| distance < b.1) {
                        best = Some((c, distance));
                    }
                }
            }
            if let Some((c, _)) = best {
                mates.push(c);
            }
        }
        Ok(mates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blosum::AnchorBlosumMatch;
    use crate::genes::VjAnchorTemplate;
    use crate::read::test_util::test_read;
    use crate::read::{cigar_from_str, FLAG_FIRST, FLAG_MATE_REVERSE, FLAG_PAIRED, FLAG_UNMAPPED};
    use std::sync::Arc;

    const MAX_FRAGMENT_LENGTH: i32 = 1000;

    struct MockSearcher {
        m: Option<AnchorBlosumMatch>,
    }
    impl AnchorSearcher for MockSearcher {
        fn search(&self, _: &[u8], _: &[VjGeneType], _: SearchMode, _: i32, _: i32) -> Result<Option<AnchorBlosumMatch>> {
            Ok(self.m.clone())
        }
    }

    fn mapped_read(start: i32, cigar: &str, seq: &str) -> CiderRead {
        let mut r = test_read();
        r.tid = 0;
        r.contig = Arc::from("X");
        r.alignment_start = start;
        r.cigar = cigar_from_str(cigar);
        r.bases = seq.as_bytes().to_vec();
        r.quals = vec![37; seq.len()];
        r.flags = FLAG_PAIRED;
        r
    }

    // CiderReadScreenerTest.testExtrapolateReadOffsetAtRefPosition
    #[test]
    fn extrapolate_read_offset() {
        let r = mapped_read(1001, "30S50M20S", "GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT");
        for java in [false, true] {
            compat::with_java(java, || {
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1001), 30);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1050), 79);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1001 - 30), 0);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1050 + 20), 99);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1050 + 30), 109);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 1001 - 32), -2);
            });
        }
    }

    // CiderReadScreenerTest.testExtrapolateReadOffsetAtRefPositionSplice
    #[test]
    fn extrapolate_read_offset_splice() {
        let seq = "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCA".repeat(2);
        let r = mapped_read(100000, "20S50M7000N20M1000N30M80S", &seq);
        for java in [false, true] {
            compat::with_java(java, || {
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 100000 - 30), -10);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 100000 - 20), 0);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 100020), 40);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 100150), 170);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 100250), 270);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 107000), 20);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 107060), 80);
                assert_eq!(extrapolate_read_offset_at_ref_position(&r, 107100), 120);
            });
        }
    }

    #[test]
    fn extrapolate_read_offset_in_deletion() {
        // CI4: 60M10D90M; a position 3 bases into the deletion is closest to the first block
        let r = mapped_read(1000, "60M10D90M", &"A".repeat(150));
        assert_eq!(compat::with_java(false, || extrapolate_read_offset_at_ref_position(&r, 1062)), 62);
        assert_eq!(compat::with_java(true, || extrapolate_read_offset_at_ref_position(&r, 1062)), 52);
    }

    // CiderReadScreenerTest.testFindAnchorPositionRNA
    #[test]
    fn find_anchor_position_rna() {
        let mut r = mapped_read(
            106322274,
            "49M7085N20M1389N16M66S",
            "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCAAGCCATTACTATAAAATTTCGGTCTCGCACAGTAATACACAGCCGTGTCTT",
        );
        r.contig = Arc::from("14");
        r.name = "read1".into();
        r.flags |= FLAG_FIRST;
        let r: ReadRef = Arc::new(r);
        let template = VjAnchorTemplate::new(
            VjGeneType::Ighj,
            "IGHJ1",
            "01",
            Some(GenomicLocation::new("14", 106330701, 106330840, Strand::Reverse)),
            "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
            Some(GenomicLocation::new("14", 106330801, 106330830, Strand::Reverse)),
        );
        let ds = GeneDatastore::new(vec![template], Vec::new()).unwrap();
        let searcher = MockSearcher { m: None };
        let screener = CiderReadScreener::new(&ds, &searcher, 6, MAX_FRAGMENT_LENGTH);
        let mapped = mapped_region(&r);
        let loc = VjAnchorGenomeLocation { gene_type: VjGeneType::Ighj, location: GenomicLocation::new("14", 106330801, 106330830, Strand::Reverse) };
        let c = screener.try_match_to_anchor_location(&r, &mapped, &loc, 0).unwrap().unwrap();
        assert_eq!(c.anchor_offset_start, 68);
        assert_eq!(c.anchor_offset_end, 98);
    }

    fn unmapped_with_mate() -> CiderRead {
        let mut r = test_read();
        r.mate_tid = 0;
        r.mate_contig = Arc::from("1");
        r.bases = b"GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT".to_vec();
        r.flags = FLAG_UNMAPPED | FLAG_PAIRED;
        r
    }

    fn set_mate(r: &mut CiderRead, start: i32, reverse: bool) {
        r.mate_alignment_start = start;
        if reverse { r.flags |= FLAG_MATE_REVERSE } else { r.flags &= !FLAG_MATE_REVERSE }
    }

    // CiderReadScreenerTest.isUnamppedReadRelevantToAnchorLocV
    #[test]
    fn unmapped_read_relevant_to_anchor_v() {
        let (s, e) = (10000, 10030);
        let mut r = unmapped_with_mate();
        let loc = VjAnchorGenomeLocation { gene_type: VjGeneType::Trav, location: GenomicLocation::new("1", s, e, Strand::Forward) };
        set_mate(&mut r, s - 500, false);
        assert!(is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        let rl = r.read_length();
        set_mate(&mut r, s - MAX_FRAGMENT_LENGTH - rl, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));

        let loc = VjAnchorGenomeLocation { gene_type: VjGeneType::Trav, location: GenomicLocation::new("1", s, e, Strand::Reverse) };
        set_mate(&mut r, e + 500, true);
        assert!(is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + MAX_FRAGMENT_LENGTH, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
    }

    // CiderReadScreenerTest.isUnamppedReadRelevantToAnchorLocJ
    #[test]
    fn unmapped_read_relevant_to_anchor_j() {
        let (s, e) = (10000, 10030);
        let mut r = unmapped_with_mate();
        let loc = VjAnchorGenomeLocation { gene_type: VjGeneType::Traj, location: GenomicLocation::new("1", s, e, Strand::Forward) };
        set_mate(&mut r, e + 500, true);
        assert!(is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + MAX_FRAGMENT_LENGTH, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));

        let loc = VjAnchorGenomeLocation { gene_type: VjGeneType::Traj, location: GenomicLocation::new("1", s, e, Strand::Reverse) };
        set_mate(&mut r, s - 500, false);
        assert!(is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        let rl = r.read_length();
        set_mate(&mut r, s + MAX_FRAGMENT_LENGTH - rl, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(!is_unmapped_read_relevant_to_anchor_loc(&r, &loc, MAX_FRAGMENT_LENGTH));
    }

    // CiderReadScreenerTest.isUnamppedReadRelevantToConstantRegion
    #[test]
    fn unmapped_read_relevant_to_constant_region() {
        let (s, e) = (10000, 10030);
        let mut r = unmapped_with_mate();
        let region = GenomicLocation::new("1", s, e, Strand::Forward);
        set_mate(&mut r, e + 500, true);
        assert!(is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + MAX_FRAGMENT_LENGTH, true);
        assert!(!is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(!is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));

        let region = GenomicLocation::new("1", s, e, Strand::Reverse);
        set_mate(&mut r, s - 500, false);
        assert!(is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, e + 500, false);
        assert!(is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        let rl = r.read_length();
        set_mate(&mut r, s - MAX_FRAGMENT_LENGTH - rl, false);
        assert!(!is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
        set_mate(&mut r, s - 500, true);
        assert!(!is_unmapped_read_relevant_to_constant_region(&r, &region, MAX_FRAGMENT_LENGTH));
    }

    #[test]
    fn boundary_positions() {
        // read start backs out a left soft clip; read end adds a right soft clip and deletions
        assert_eq!(read_boundary_position(100, "10S90M", true), 90);
        assert_eq!(read_boundary_position(100, "90M10S", true), 100);
        assert_eq!(read_boundary_position(100, "10S90M", false), 189);
        assert_eq!(read_boundary_position(100, "80M5D10M10S", false), 204);
        assert_eq!(read_boundary_position(100, "5H10S90M", true), 90);
    }

    #[test]
    fn adapter_trim() {
        // forward read at 1000 (100M), its reverse mate's 5' end (unclipped end) at 1089: 10 adapter bases
        let mut r = mapped_read(1000, "100M", &"A".repeat(100));
        r.mate_tid = 0;
        r.mate_contig = r.contig.clone();
        r.mate_alignment_start = 1000;
        r.flags |= FLAG_MATE_REVERSE;
        r.mate_cigar = Some("90M".into());
        for java in [false, true] {
            assert_eq!(compat::with_java(java, || adapter_sequence_trim(&r)), (0, 10));
        }
        r.mate_cigar = None;
        assert_eq!(adapter_sequence_trim(&r), (0, 0));
    }
}
