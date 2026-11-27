//! Structural variants: StructuralVariantFactory (breakends paired by MATEID, single breakends), the legs,
//! and VariantContextCollection (the SV VCF records kept sorted by contig, position and ID, as the
//! TreeSet with htsjdk's VariantContextComparator keeps them).

use std::collections::HashMap;

use hmf_common::variant_context::{Value, VcfHeader, VcfRecord};

use crate::types::*;

#[derive(Clone, Debug, PartialEq)]
pub struct SvLeg {
    pub chromosome: String,
    pub position: i32,
    pub orientation: i8,
    pub start_offset: i32,
    pub end_offset: i32,
    pub homology: String,
    pub allele_frequency: Option<f64>,
    pub inexact_homology_offset_start: Option<i32>,
    pub inexact_homology_offset_end: Option<i32>,
    pub tumor_variant_fragment_count: Option<i32>,
    pub tumor_reference_fragment_count: Option<i32>,
    pub normal_variant_fragment_count: Option<i32>,
    pub normal_reference_fragment_count: Option<i32>,
    pub anchoring_support_distance: i32,
}

impl SvLeg {
    pub fn cna_position(&self) -> i32 {
        if self.orientation == -1 { self.position } else { self.position + 1 }
    }
}

#[derive(Clone, Debug)]
pub struct StructuralVariant {
    pub id: String,
    pub mate_id: Option<String>,
    pub start: SvLeg,
    pub end: Option<SvLeg>,
    pub sv_type: SvType,
    pub filter: String,
    pub qual: f64,
    pub hotspot: bool,
    pub insert_sequence: String,
    /// indices of the records in the collection the variant was built from
    pub start_context: usize,
    pub end_context: Option<usize>,
}

impl StructuralVariant {
    pub fn is_filtered(&self) -> bool {
        self.filter != "PASS"
    }
    pub fn chromosome(&self, start: bool) -> Option<&str> {
        if start { Some(&self.start.chromosome) } else { self.end.as_ref().map(|e| e.chromosome.as_str()) }
    }
    pub fn position(&self, start: bool) -> Option<i32> {
        if start { Some(self.start.position) } else { self.end.as_ref().map(|e| e.position) }
    }
}

/// VariantAltInsertCoords.BREAKEND_REGEX "^(.*)([\[\]])(.+)[\[\]](.*)$": (group 1, group 2, group 3, group 4)
fn breakend_groups(alt: &str) -> Option<(&str, char, &str, &str)> {
    let b = alt.as_bytes();
    let last = b.iter().rposition(|&c| c == b'[' || c == b']')?;
    if last == 0 {
        return None;
    }
    // group 1 is greedy: the last bracket before `last - 1` that leaves a non-empty group 3
    let first = b[..last - 1].iter().rposition(|&c| c == b'[' || c == b']')?;
    Some((&alt[..first], b[first] as char, &alt[first + 1..last], &alt[last + 1..]))
}

pub fn is_single_breakend(alt: &str) -> bool {
    let b = alt.as_bytes();
    !b.is_empty() && (b[0] == b'.' || b[b.len() - 1] == b'.')
}

/// SvUtils.formSvType
fn form_sv_type(chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: i8, orient_end: i8, has_inserted: bool) -> SvType {
    if chr_start != chr_end {
        return SvType::Bnd;
    }
    if orient_start != orient_end {
        let diff = (pos_start - pos_end).abs();
        if diff == 1 && has_inserted {
            return SvType::Ins;
        }
        if diff == 0 {
            return SvType::Dup;
        }
        let first_is_lower = pos_start < pos_end;
        return if first_is_lower == (orient_start == 1) { SvType::Del } else { SvType::Dup };
    }
    SvType::Inv
}

/// Which genotypes give the tumor and normal fragment counts
#[derive(Clone, Copy, Debug)]
pub struct GenotypeOrdinals {
    pub reference: i32,
    pub tumor: i32,
}

pub struct SvFactory<'a> {
    header: &'a VcfHeader,
    ordinals: Option<GenotypeOrdinals>,
    unmatched: HashMap<String, usize>,
    /// IDs in the order they were added to `unmatched` (Java iterates a HashMap: String hash order)
    pub results: Vec<StructuralVariant>,
}

fn int_attr(g: &hmf_common::variant_context::Genotype, key: &str) -> i32 {
    g.attr_int(key, 0)
}

impl<'a> SvFactory<'a> {
    pub fn new(header: &'a VcfHeader, ordinals: Option<GenotypeOrdinals>) -> SvFactory<'a> {
        SvFactory { header, ordinals, unmatched: HashMap::new(), results: Vec::new() }
    }

    /// StructuralVariantFactory.addVariantContext; `records` holds every record, `index` the one added
    pub fn add(&mut self, records: &mut [VcfRecord], index: usize) {
        let alt = records[index].alt().to_string();
        if is_single_breakend(&alt) {
            let sv = self.create_single_breakend(records, index);
            self.results.push(sv);
        } else {
            let mate = records[index].info.get("MATEID").map(|v| v.java_string());
            if let Some(m) = mate.as_ref().and_then(|m| self.unmatched.remove(m)) {
                let sv = self.create_sv(records, m, index);
                self.results.push(sv);
            } else {
                self.unmatched.insert(records[index].id.clone(), index);
            }
        }
    }

    /// the records left without a mate
    pub fn unmatched(&self) -> Vec<usize> {
        let mut keys: Vec<(&String, &usize)> = self.unmatched.iter().collect();
        keys.sort_by(|a, b| a.0.cmp(b.0));
        keys.into_iter().map(|(_, i)| *i).collect()
    }

    fn leg(&self, records: &mut [VcfRecord], index: usize, ignore_refpair: bool, orientation: i8) -> SvLeg {
        let rec = &mut records[index];
        let (mut ci_left, mut ci_right) = (0, 0);
        if rec.info.has("CIPOS") {
            let cipos = rec.info.as_int_list("CIPOS", 0);
            if cipos.len() == 2 {
                ci_left = cipos[0];
                ci_right = cipos[1];
            }
        }
        let anchor = rec.info.as_int_list("SEGALEN", 0).into_iter().max().unwrap_or(0);
        let (mut ih_start, mut ih_end) = (Some(0), Some(0));
        if rec.info.has("IHOMPOS") {
            let ihom = rec.info.as_int_list("IHOMPOS", 0);
            if ihom.len() == 2 {
                ih_start = Some(ihom[0]);
                ih_end = Some(ihom[1]);
            }
        }
        let genotypes = rec.genotypes(self.header);
        let (reference, tumor) = match self.ordinals {
            Some(o) => (o.reference, o.tumor),
            None => {
                if genotypes.len() == 1 {
                    (-1, 0)
                } else {
                    (0, 1)
                }
            }
        };
        let mut leg = SvLeg {
            chromosome: rec.chrom.clone(),
            position: rec.pos,
            orientation,
            start_offset: ci_left,
            end_offset: ci_right,
            homology: String::new(),
            allele_frequency: None,
            inexact_homology_offset_start: ih_start,
            inexact_homology_offset_end: ih_end,
            tumor_variant_fragment_count: None,
            tumor_reference_fragment_count: None,
            normal_variant_fragment_count: None,
            normal_reference_fragment_count: None,
            anchoring_support_distance: anchor,
        };
        let genotypes = match &rec.genotypes {
            hmf_common::variant_context::Genotypes::Decoded(g) => g,
            _ => unreachable!(),
        };
        if reference >= 0 {
            if let Some(g) = genotypes.get(reference as usize) {
                let total = int_attr(g, "VF");
                let refs = int_attr(g, "REF");
                let pairs = int_attr(g, "REFPAIR");
                leg.normal_variant_fragment_count = Some(total);
                leg.normal_reference_fragment_count = Some(refs + if ignore_refpair { 0 } else { pairs });
            }
        }
        if tumor >= 0 {
            if let Some(g) = genotypes.get(tumor as usize) {
                let total = int_attr(g, "VF");
                let refs = int_attr(g, "REF");
                let pairs = int_attr(g, "REFPAIR");
                leg.tumor_variant_fragment_count = Some(total);
                leg.tumor_reference_fragment_count = Some(refs + if ignore_refpair { 0 } else { pairs });
                leg.allele_frequency = Some(g.attr_double("AF", 0.0));
            }
        }
        leg
    }

    /// StructuralVariantFactory.filters
    fn filters(a: &VcfRecord, b: Option<&VcfRecord>) -> String {
        let mut set: Vec<String> = Vec::new();
        for r in std::iter::once(a).chain(b) {
            if let Some(f) = &r.filters {
                for x in f {
                    if !set.contains(x) {
                        set.push(x.clone());
                    }
                }
            }
        }
        if set.len() > 1 {
            set.retain(|f| f != "PASS");
        }
        if set.is_empty() {
            set.push("PASS".into());
        }
        set.sort();
        set.join(";")
    }

    fn common_hotspot(rec: &VcfRecord) -> bool {
        rec.info.as_bool("HOTSPOT", false)
    }

    pub fn create_sv(&self, records: &mut [VcfRecord], start_index: usize, end_index: usize) -> StructuralVariant {
        let alt = records[start_index].alt().to_string();
        let (g1, g2, _g3, g4) = breakend_groups(&alt).unwrap_or_else(|| panic!("ALT {alt} is not in breakend notation"));
        let mut start_orient: i8 = if !g1.is_empty() { 1 } else { -1 };
        let mut end_orient: i8 = if g2 == ']' { 1 } else { -1 };
        let inserted = if !g1.is_empty() { g1[1..].to_string() } else { g4[..g4.len() - 1].to_string() };
        let sv_type_attr = records[start_index].info.get("SVTYPE").map(|v| v.java_string());
        let (mut ls, mut le) = (start_index, end_index);
        if records[start_index].chrom == records[end_index].chrom
            && start_orient != end_orient
            && records[start_index].pos == records[end_index].pos
            && start_orient == 1
        {
            ls = end_index;
            le = start_index;
            start_orient = -1;
            end_orient = 1;
        }
        let sv_type = match sv_type_attr.as_deref().and_then(SvType::from_name) {
            Some(t) => t,
            None => form_sv_type(
                &records[ls].chrom,
                &records[le].chrom,
                records[ls].pos,
                records[le].pos,
                start_orient,
                end_orient,
                !inserted.is_empty(),
            ),
        };
        let indel_length = if matches!(sv_type, SvType::Del | SvType::Dup | SvType::Ins) {
            (records[start_index].pos - records[end_index].pos).abs()
        } else {
            0
        };
        let small = matches!(sv_type, SvType::Del | SvType::Dup | SvType::Ins) && indel_length <= 1000;
        let mut start_leg = self.leg(records, ls, small, start_orient);
        start_leg.homology = records[ls].info.as_string("HOMSEQ", "");
        let mut end_leg = self.leg(records, le, small, end_orient);
        end_leg.homology = records[le].info.as_string("HOMSEQ", "");
        let src = &records[start_index];
        StructuralVariant {
            id: src.id.clone(),
            mate_id: Some(records[le].id.clone()),
            start: start_leg,
            end: Some(end_leg),
            sv_type,
            filter: Self::filters(&records[ls], Some(&records[le])),
            qual: src.phred_scaled_qual(),
            hotspot: Self::common_hotspot(src),
            insert_sequence: inserted,
            start_context: ls,
            end_context: Some(le),
        }
    }

    pub fn create_single_breakend(&self, records: &mut [VcfRecord], index: usize) -> StructuralVariant {
        let alt = records[index].alt().to_string();
        let orientation: i8 = if alt.starts_with('.') { -1 } else { 1 };
        let ref_len = records[index].ref_allele.len();
        let inserted = if orientation == -1 { alt[1..alt.len() - ref_len].to_string() } else { alt[ref_len..alt.len() - 1].to_string() };
        let leg = self.leg(records, index, false, orientation);
        let rec = &records[index];
        StructuralVariant {
            id: rec.id.clone(),
            mate_id: None,
            start: leg,
            end: None,
            sv_type: if rec.info.has("INFERRED") { SvType::Inf } else { SvType::Sgl },
            filter: Self::filters(rec, None),
            qual: rec.phred_scaled_qual(),
            hotspot: Self::common_hotspot(rec),
            insert_sequence: inserted,
            start_context: index,
            end_context: None,
        }
    }
}

/// VariantContextCollection: records sorted by (contig index, position, ID) without duplicates
pub struct SvCollection {
    pub header: VcfHeader,
    contig_index: HashMap<String, usize>,
    pub records: Vec<VcfRecord>,
}

impl SvCollection {
    pub fn new(header: VcfHeader) -> SvCollection {
        let contig_index = header.contigs().into_iter().enumerate().map(|(i, c)| (c, i)).collect();
        SvCollection { header, contig_index, records: Vec::new() }
    }

    fn key(&self, r: &VcfRecord) -> (usize, i32) {
        let c = *self.contig_index.get(&r.chrom).unwrap_or_else(|| panic!("contig {} not in the VCF header's dictionary", r.chrom));
        (c, r.pos)
    }

    /// TreeSet.add: a record equal in contig, position and ID to one present is not added
    pub fn add(&mut self, rec: VcfRecord) {
        let k = self.key(&rec);
        let pos = self.records.partition_point(|x| {
            let xk = self.key(x);
            xk < k || (xk == k && hmf_common::java_text::compare(&x.id, &rec.id) == std::cmp::Ordering::Less)
        });
        if pos < self.records.len() && self.key(&self.records[pos]) == k && self.records[pos].id == rec.id {
            return;
        }
        self.records.insert(pos, rec);
    }

    /// VariantContextCollection.variants: the SVs of the passing and inferred records on human chromosomes
    pub fn variants(&mut self) -> Vec<StructuralVariant> {
        let header = self.header.clone();
        let mut factory = SvFactory::new(&header, None);
        for i in 0..self.records.len() {
            let r = &self.records[i];
            let pass = !r.is_filtered() || r.info.has("INFERRED");
            if pass && Chromosome::contains(&r.chrom) {
                factory.add(&mut self.records, i);
            }
        }
        factory.results
    }
}

/// a breakend added by PURPLE where a copy number change has no SV (SomaticSvCache.infer)
pub fn inferred_record(chromosome: &str, start: i32, min_start: i32, max_start: i32, increasing: bool, id: usize) -> VcfRecord {
    let (alt, position) = if increasing { (".N", start) } else { ("N.", start - 1) };
    let lower = (-500i64).min(min_start as i64 - start as i64);
    let upper = 500i64.max(max_start as i64 - start as i64);
    let mut info = hmf_common::variant_context::Attributes::default();
    info.put("CIPOS", Value::IntList(vec![lower as i32, upper as i32]));
    info.put("SVTYPE", Value::Str("BND".into()));
    info.put("INFERRED", Value::Bool(true));
    VcfRecord {
        chrom: chromosome.to_string(),
        pos: position,
        id: format!("purple_{id}"),
        ref_allele: "N".into(),
        alts: vec![alt.into()],
        log10_p_error: hmf_common::variant_context::NO_LOG10_PERROR,
        filters: Some(vec!["INFERRED".into()]),
        info,
        genotypes: hmf_common::variant_context::Genotypes::Decoded(Vec::new()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breakend_notation() {
        assert_eq!(breakend_groups("ACGT[chr2:123["), Some(("ACGT", '[', "chr2:123", "")));
        assert_eq!(breakend_groups("]chr2:123]ACGT"), Some(("", ']', "chr2:123", "ACGT")));
        assert_eq!(breakend_groups("G[chr1:5023347["), Some(("G", '[', "chr1:5023347", "")));
        assert!(breakend_groups("ACGT").is_none());
        assert!(is_single_breakend(".ACGT") && is_single_breakend("ACGT.") && !is_single_breakend("A[1:2["));
    }

    #[test]
    fn sv_types() {
        assert_eq!(form_sv_type("1", "2", 1, 2, 1, 1, false), SvType::Bnd);
        assert_eq!(form_sv_type("1", "1", 100, 200, 1, -1, false), SvType::Del);
        assert_eq!(form_sv_type("1", "1", 100, 200, -1, 1, false), SvType::Dup);
        assert_eq!(form_sv_type("1", "1", 100, 101, 1, -1, true), SvType::Ins);
        assert_eq!(form_sv_type("1", "1", 100, 200, 1, 1, false), SvType::Inv);
    }
}
