//! The tumor side with PURPLE: SomaticVariantAnnotation (the somatic variants in the HLA genes' coding
//! regions, and the alleles their fragments support), the variant fragments of BamRecordReader, and
//! CopyNumberAssignment.

use anyhow::{bail, Context, Result};

use crate::bam::{filter_low_qual_read, Bam, Read, ReadCtx};
use crate::config::LilacConfig;
use crate::constants::SPLICE_VARIANT_BUFFER;
use crate::coverage::AlleleCoverage;
use crate::fragment::{merge_fragments_by_id, Factory, Fragment};
use crate::hla::{AId, HlaGene};
use crate::indel::Indel;
use crate::qc::CodingEffect;
use crate::reader::GeneReader;
use crate::refdata::{nucleotide_locus, RefData, Transcript};
use crate::seq::SeqLoci;
use crate::tok::{self, WILD};

#[derive(Clone, Debug)]
pub struct SomaticVariant {
    pub gene: Option<HlaGene>,
    pub chromosome: String,
    pub position: i32,
    pub ref_: String,
    pub alt: String,
    pub filter: String,
    pub effect: CodingEffect,
}

impl SomaticVariant {
    pub fn display(&self) -> String {
        format!(
            "gene({}) pos({}:{}) variant({}>{} {})",
            self.gene.map(|g| g.long_name()).unwrap_or_else(|| "null".to_string()),
            self.chromosome,
            self.position,
            self.ref_,
            self.alt,
            self.effect.name()
        )
    }
}

/// The VCF records LILAC reads: PASS (or unfiltered), with the PAVE IMPACT annotation
fn read_vcf(path: &str, hla_chr: &str) -> Result<Vec<(String, i32, String, String, String, Vec<String>)>> {
    if !std::path::Path::new(path).exists() {
        bail!("java.lang.NullPointerException (somatic VCF {path} not found)");
    }
    let mut out = Vec::new();
    let prefix = format!("{hla_chr}\t");
    hmf_common::java_text::for_each_line(path, |line| {
        if line.starts_with('#') || !line.starts_with(&prefix) {
            return Ok(());
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Ok(());
        }
        let filter = f[6];
        // isFiltered() && !getFilters().contains(PASS): kept when unfiltered (".") or PASS
        let filters: Vec<&str> = if filter == "." { Vec::new() } else { filter.split(';').collect() };
        if !filters.is_empty() && !filters.contains(&"PASS") {
            return Ok(());
        }
        let pos = hmf_common::java_text::parse_int(f[1]).context("invalid VCF position")?;
        let mut impact: Vec<String> = Vec::new();
        for kv in f[7].split(';') {
            if let Some(v) = kv.strip_prefix("IMPACT=") {
                impact = v.split(',').map(|x| x.to_string()).collect();
            }
        }
        let display_filter = if filters.is_empty() || filters == ["PASS"] { "PASS".to_string() } else { filters.join(";") };
        out.push((f[0].to_string(), pos, f[3].to_string(), f[4].to_string(), display_filter, impact));
        Ok(())
    })?;
    Ok(out)
}

pub struct VariantAnnotation {
    pub variants: Vec<SomaticVariant>,
    gene_loci: Vec<(Option<HlaGene>, Vec<i32>)>,
}

impl VariantAnnotation {
    /// new SomaticVariantAnnotation(config, transcripts)
    pub fn new(config: &LilacConfig, rd: &RefData) -> Result<VariantAnnotation> {
        let mut va = VariantAnnotation { variants: Vec::new(), gene_loci: Vec::new() };
        if config.somatic_variants_file.is_empty() {
            return Ok(va);
        }
        let hla_chr = &config.params.hla_chr;
        let transcripts: Vec<&Transcript> = rd.transcripts.iter().map(|x| &x.1).collect();
        let mut loaded = Vec::new();
        for (chr, pos, r, a, filter, impact) in read_vcf(&config.somatic_variants_file, hla_chr)? {
            let (gene_name, effect) = if impact.len() == 10 {
                (impact[0].clone(), CodingEffect::parse(&impact[3]).context("invalid coding effect")?)
            } else {
                (String::new(), CodingEffect::Undefined)
            };
            // inHlaCodingRegion
            let end = pos + r.len() as i32 - 1;
            let mut in_region = false;
            for t in &transcripts {
                if !(pos <= t.coding_end && end >= t.coding_start) {
                    continue;
                }
                if t.exons.iter().any(|e| pos <= e.end && end >= e.start) || effect == CodingEffect::Splice {
                    in_region = true;
                    break;
                }
            }
            if !in_region {
                continue;
            }
            loaded.push(SomaticVariant { gene: HlaGene::from_string(&gene_name), chromosome: chr, position: pos, ref_: r, alt: a, filter, effect });
        }
        crate::info!("loaded {} HLA variants from file: {}", loaded.len(), config.somatic_variants_file);
        for v in loaded {
            let locus = nucleotide_locus(transcripts.iter().copied(), v.position);
            if locus < 0 {
                continue;
            }
            match va.gene_loci.iter_mut().find(|x| x.0 == v.gene) {
                Some(x) => x.1.push(locus / 3),
                None => va.gene_loci.push((v.gene, vec![locus / 3])),
            }
            va.variants.push(v);
        }
        Ok(va)
    }

    /// assignAlleleCoverage: the winning alleles most consistent with the variant's fragments
    pub fn assign_allele_coverage(&self, rd: &RefData, config: &LilacConfig, variant: &SomaticVariant, frags: &mut [Fragment], winners: &[&SeqLoci]) -> Vec<AlleleCoverage> {
        let mut out: Vec<AlleleCoverage> = Vec::new();
        if frags.is_empty() {
            return out;
        }
        for f in frags.iter_mut() {
            f.remove_low_qual_bases();
        }
        for f in frags.iter_mut() {
            f.build_amino_acids();
        }
        let loci = self.gene_loci.iter().find(|x| x.0 == variant.gene).map(|x| x.1.as_slice());
        for s in winners {
            let Some(vl) = loci.filter(|v| !v.is_empty()) else {
                crate::error!("no sequences found for variant({})", variant.display());
                continue;
            };
            let mut support = 0;
            for f in frags.iter() {
                let mut ok = true;
                let mut n = 0;
                let (lo, hi) = (f.min_aa_locus(), f.max_aa_locus());
                let mut locus = lo;
                while locus <= hi {
                    if locus as usize >= s.len() {
                        break;
                    }
                    if vl.contains(&locus) {
                        locus += 1;
                        continue;
                    }
                    let aa = match f.aas.get(locus) {
                        Some(a) => a,
                        None => {
                            let l = f.low_qual_aa(locus, false, &config.params);
                            if l == tok::EMPTY {
                                locus += 1;
                                continue;
                            }
                            l
                        }
                    };
                    let sq = s.at(locus as usize);
                    if tok::to_string(aa) != tok::to_string(sq) && sq != WILD {
                        ok = false;
                        break;
                    }
                    n += 1;
                    locus += 1;
                }
                if ok && n > 0 {
                    support += 1;
                }
            }
            if support > 0 {
                if s.has_wildcards() {
                    out.push(AlleleCoverage::new(s.allele, 0, 0.0, support as f64));
                } else {
                    out.push(AlleleCoverage::new(s.allele, support, 0.0, 0.0));
                }
            }
        }
        let _ = rd;
        let max_unique = out.iter().map(|x| x.unique).max().unwrap_or(0);
        if max_unique > 0 && max_unique as usize >= frags.len() / 2 {
            return out.into_iter().filter(|x| x.unique == max_unique).collect();
        }
        let max_wild = out.iter().map(|x| x.wild as i32).max().unwrap_or(0);
        if max_wild > 0 {
            return out.into_iter().filter(|x| x.wild as i32 == max_wild).collect();
        }
        Vec::new()
    }
}

/// regionCoversVariant
fn region_covers(region: (i32, i32), v: &SomaticVariant) -> bool {
    if v.effect == CodingEffect::Splice {
        return (v.position - region.0).abs() <= SPLICE_VARIANT_BUFFER || (v.position - region.1).abs() <= SPLICE_VARIANT_BUFFER;
    }
    v.position >= region.0 && v.position <= region.1
}

/// recordContainsVariant
fn record_contains_variant(v: &SomaticVariant, r: &Read) -> bool {
    if v.alt.len() != v.ref_.len() {
        let expected = Indel::of(&v.chromosome, v.position, &v.ref_, &v.alt);
        return r.valid_indels.iter().any(|x| x.matches(&expected));
    }
    for (i, b) in v.alt.bytes().enumerate() {
        let pos = v.position + i as i32;
        let idx = r.record.read_pos_at_ref_pos(pos, false) - 1;
        if idx < 0 {
            return false;
        }
        if r.record.seq.get(idx as usize) != Some(&b) {
            return false;
        }
    }
    true
}

/// BamRecordReader.findVariantFragments
pub fn variant_fragments(reader: &mut GeneReader, bam: &mut Bam, factory: &Factory, ctx: &ReadCtx, v: &SomaticVariant) -> Result<Vec<Fragment>> {
    let regions: Vec<(HlaGene, i8, Vec<(i32, i32)>)> = reader.regions.iter().map(|r| (r.gene, r.strand, r.coding_regions.clone())).collect();
    for (gene, strand, coding) in regions {
        for region in coding.iter().copied() {
            if !region_covers(region, v) {
                continue;
            }
            let slice = if v.effect == CodingEffect::Splice { (region.0 - SPLICE_VARIANT_BUFFER, region.1 + SPLICE_VARIANT_BUFFER) } else { region };
            // findVariantRecords: the reads at the variant position
            let mut reads: Vec<Read> = Vec::new();
            for rec in bam.slice(&v.chromosome, v.position, v.position)? {
                if filter_low_qual_read(&rec, ctx.params) {
                    continue;
                }
                if reader.both_ends_in_range(&rec, &ctx.params.hla_chr) {
                    reads.push(Read::create(ctx, slice, &rec, true, true)?);
                } else {
                    reader.filtered_count += 1;
                }
            }
            let mut coding_reads: Vec<Read> = Vec::new();
            for r in reads {
                if !record_contains_variant(v, &r) {
                    continue;
                }
                // a record equal to one taken already (SAMRecord.hashCode covers its content)
                if coding_reads.iter().any(|x| x.record.name == r.record.name && x.record.flags == r.record.flags && x.record.start == r.record.start && x.record.cigar == r.record.cigar && x.record.seq == r.record.seq) {
                    continue;
                }
                coding_reads.push(r);
            }
            let mut frags: Vec<Fragment> = Vec::new();
            for r in &coding_reads {
                if let Some(f) = factory.create_alignment_fragments(r, gene, strand)? {
                    frags.push(f);
                }
            }
            // queryMateFragments
            for r in &coding_reads {
                let Some(mate) = bam.query_mate(&r.record)? else { continue };
                let Some(cr) = coding.iter().copied().find(|c| c.0 <= mate.end && c.1 >= mate.start) else { continue };
                let mr = Read::create(ctx, cr, &mate, true, true)?;
                if let Some(f) = factory.create_alignment_fragments(&mr, gene, strand)? {
                    frags.push(f);
                }
            }
            return Ok(merge_fragments_by_id(frags));
        }
    }
    Ok(Vec::new())
}

/// The PURPLE gene copy numbers of the HLA genes: (gene, minCopyNumber, minMinorAlleleCopyNumber)
pub fn load_copy_numbers(path: &str, genes: &[HlaGene]) -> Vec<(HlaGene, f64, f64)> {
    let lines = match hmf_common::java_text::read_lines(path) {
        Ok(l) => l,
        Err(e) => {
            crate::error!("failed to read gene copy number file({path}): {e}");
            return Vec::new();
        }
    };
    if lines.is_empty() {
        return Vec::new();
    }
    let idx = hmf_common::java_text::fields_index_map(&lines[0], '\t');
    let (Some(gi), Some(mi), Some(ni)) = (idx.get("gene"), idx.get("minCopyNumber"), idx.get("minMinorAlleleCopyNumber")) else {
        crate::error!("failed to read gene copy number file({path}): missing columns");
        return Vec::new();
    };
    let mut out = Vec::new();
    for l in &lines[1..] {
        let f = hmf_common::java_text::split(l, '\t');
        let Some(g) = f.get(*gi).and_then(|x| HlaGene::from_string(x)) else { continue };
        if !genes.contains(&g) {
            continue;
        }
        let min = f.get(*mi).and_then(|x| hmf_common::java_text::parse_double(x)).unwrap_or(0.0);
        let minor = f.get(*ni).and_then(|x| hmf_common::java_text::parse_double(x)).unwrap_or(0.0);
        out.push((g, min, minor));
    }
    out
}

/// CopyNumberAssignment.assign: each allele's tumor copy number from its gene's major and minor copy numbers
pub fn assign_copy_numbers(rd: &RefData, cn: &[(HlaGene, f64, f64)], reference: &[AlleleCoverage], tumor: &[AlleleCoverage], out: &mut Vec<f64>) {
    crate::info!("calculating tumor copy number of winning alleles");
    let n = rd.gene_cache.expect_allele_count;
    if reference.len() != n || tumor.len() != n {
        return;
    }
    if cn.is_empty() {
        return;
    }
    out.clear();
    let mut i = 0;
    while i < n {
        let (r1, r2, t1, t2) = (&reference[i], &reference[i + 1], &tumor[i], &tumor[i + 1]);
        let gene = rd.alleles.gene(r1.allele);
        match cn.iter().find(|x| x.0 == gene) {
            None => {
                crate::warn!("missing gene({}) copy number data", gene.long_name());
                out.push(0.0);
                out.push(0.0);
            }
            Some(c) => {
                let minor = c.2;
                let major = c.1 - minor;
                if rd.alleles.name(r1.allele) == rd.alleles.name(r2.allele) {
                    out.push(major);
                    out.push(minor);
                } else {
                    let ratio1 = t1.total / r1.total;
                    let ratio2 = t2.total / r2.total;
                    if ratio1 >= ratio2 {
                        out.push(major);
                        out.push(minor);
                    } else {
                        out.push(minor);
                        out.push(major);
                    }
                }
            }
        }
        i += 2;
    }
    let _: Option<AId> = None;
}
