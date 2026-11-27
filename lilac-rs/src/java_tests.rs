//! The Java unit tests of LILAC (hmftools/lilac/src/test), ported. Where they mock a SequenceCount or a
//! Fragment with Mockito, the port builds real ones with the same answers.

use std::sync::Arc;

use crate::bam::{Read, ReadCtx, Rec};
use crate::config::{LilacConfig, Params};
use crate::coverage::{self, AlleleCoverage, ComplexBuilder, ComplexCoverage, FragAlleles, Mapper};
use crate::evidence::{self, PhasedEvidence};
use crate::fragment::{self, Fragment, Nuc, Scope};
use crate::hla::{AId, Alleles, GeneSelector, HlaGene};
use crate::qc::{self, AminoAcidQC, CodingEffect, Haplotype, SomaticCodingCount};
use crate::refdata::{self, GeneCache, GeneEnrichment, RefData};
use crate::seq::{self, MatchType, SeqLoci, SequenceCount};
use crate::tok::{self, Tok};

fn t(s: &str) -> Tok {
    tok::intern(s)
}
fn toks(v: &[&str]) -> Vec<Tok> {
    v.iter().map(|x| t(x)).collect()
}

/// SamRecordTestUtils.createSamRecord: paired, proper pair, first in pair, base qual 37
fn sam_record(id: &str, chr: &str, start: i32, bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, reversed: bool) -> Rec {
    let mut ops: Vec<(u8, u32)> = Vec::new();
    let mut n = 0u32;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            ops.push((c, n));
            n = 0;
        }
    }
    let ref_len: i32 = ops.iter().filter(|c| matches!(c.0, b'M' | b'D' | b'N' | b'=' | b'X')).map(|c| c.1 as i32).sum();
    let mut flags = 0x1 | 0x2 | 0x40;
    if reversed {
        flags |= 0x10;
    }
    Rec {
        name: id.to_string(),
        flags,
        contig: chr.to_string(),
        start,
        end: start + ref_len - 1,
        cigar: ops,
        seq: bases.as_bytes().to_vec(),
        quals: vec![37; bases.len()],
        mapq: 60,
        mate_contig: mate_chr.to_string(),
        mate_tid: 0,
        mate_start,
        isize: if chr == mate_chr { (start - mate_start).abs() as i64 } else { 0 },
        uq: None,
    }
}

fn random_bases(n: usize) -> String {
    let b = b"ACGT";
    let mut x: u64 = 12345;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b[(x >> 33) as usize % 4] as char
        })
        .collect()
}

/// LilacTestUtils.createReadRecord
fn read_record(id: &str) -> Arc<Rec> {
    Arc::new(sam_record(id, "1", 100, &random_bases(151), "151M", "1", 300, false))
}

/// Fragment.createFromQuals
fn frag_from_quals(id: &str, gene: HlaGene, genes: Vec<HlaGene>, loci: &[i32], quals: &[i8], bases: &[&str]) -> Fragment {
    let nucs: Vec<(i32, Nuc)> = loci.iter().enumerate().map(|(i, l)| (*l, Nuc { bases: t(bases[i]), low_qual: fragment::is_low_qual_byte(quals[i]), qual: quals[i] })).collect();
    Fragment::new(read_record(id), gene, genes, nucs)
}

/// LilacTestUtils.createFragment(id, gene, sequence, locusStart, locusEnd)
fn test_fragment(id: &str, gene: HlaGene, sequence: &str, start: i32, end: i32) -> Fragment {
    let loci: Vec<i32> = (start..=end).collect();
    let b = sequence.as_bytes();
    let bases: Vec<String> = if loci.len() > b.len() { Vec::new() } else { loci.iter().map(|x| (b[*x as usize] as char).to_string()).collect() };
    let refs: Vec<&str> = bases.iter().map(|x| x.as_str()).collect();
    let loci = if refs.is_empty() { Vec::new() } else { loci };
    frag_from_quals(id, gene, vec![gene], &loci, &vec![37; loci.len()], &refs)
}

fn empty_fragment(id: &str) -> Fragment {
    Fragment::new(read_record(id), HlaGene::None, Vec::new(), Vec::<(i32, Nuc)>::new())
}

fn class1_v37() -> GeneCache {
    GeneCache::new(&refdata::load_hla_transcripts(false, Some(GeneSelector::MhcClass1)))
}

fn params() -> Params {
    Params::default()
}

fn counts(loci: &[(i32, &str, u32)]) -> SequenceCount {
    let mut c = SequenceCount::new(0.0, 1);
    for (l, s, n) in loci {
        for _ in 0..*n {
            c.add(*l, t(s));
        }
    }
    c
}

/// HlaSequenceLoci.create(allele, sequence, sequence)
fn seq_loci(al: &mut Alleles, allele: &str, seq: &str) -> SeqLoci {
    let a = al.request_four_digit(allele).unwrap();
    let s = seq::create_aligned(seq, seq);
    SeqLoci::new(al, a, s)
}

fn seq_of(al: &mut Alleles, allele: &str, v: &[&str]) -> SeqLoci {
    let a = al.request_four_digit(allele).unwrap();
    SeqLoci::new(al, a, toks(v))
}

// ---------------------------------------------------------------- hla

#[test]
fn hla_allele_decode_and_reduce() {
    for s in ["A*01:01:01:01", "A*01:01:01", "A*26:22", "A*26"] {
        let p = crate::hla::parse_allele(s).unwrap();
        assert_eq!(crate::hla::allele_name(p.gene, &p.group, &p.protein, &p.synonymous, &p.synonymous_non_coding), s);
    }
    let mut al = Alleles::new();
    let eight = al.request("A*01:01:01:01").unwrap();
    assert_eq!(al.name(al.four_digit(eight)), "A*01:01");
    assert_eq!(al.name(al.group(eight)), "A*01");
    assert_eq!(refdata::four_digit_name("A*01:01:01:01").unwrap(), "A*01:01");
}

#[test]
fn hla_allele_cache() {
    let mut al = Alleles::new();
    let a1 = al.request("A*01:01:01").unwrap();
    let a2 = al.request("A*01:01:02").unwrap();
    let a3 = al.request("A*01:02:01").unwrap();
    let a4 = al.request("A*01:02:02").unwrap();
    let f1 = al.request_four_digit("A*01:01").unwrap();
    assert_eq!(f1, al.four_digit(a1));
    assert_eq!(f1, al.four_digit(a2));
    assert_eq!(al.four_digit(a3), al.four_digit(a4));
    assert_eq!(al.group(a1), al.group(a4));
    let f2 = al.request_four_digit("A*01:02").unwrap();
    assert_eq!(f2, al.four_digit(a3));
    assert_eq!(al.group(f1), al.group(a1));
}

#[test]
fn hla_gene_cache_class1() {
    let t = refdata::load_hla_transcripts(false, Some(GeneSelector::MhcClass1));
    assert_eq!(t.len(), 3);
    assert_eq!(GeneCache::new(&t).gene_names.len(), 3);
}

#[test]
fn hla_sequence_loci() {
    let mut al = Alleles::new();
    let a = al.request("A*01:01:01:01").unwrap();
    let r = "APRTSMNVEPRF";
    let s = SeqLoci::new(&mut al, a, seq::create_aligned(r, r));
    assert_eq!(s.len(), 12);
    assert_eq!(s.sequence_str(), r);
    assert!(!s.has_indels() && !s.has_inserts && !s.has_deletes);

    let p = al.request_four_digit("A*01:01").unwrap();
    let ins = SeqLoci::new(&mut al, p, seq::parse_reference_sequence("APRT|SET|MNVEPRF"));
    assert_eq!(ins.len(), 12);
    assert_eq!(tok::to_string(ins.at(4)), "SET");
    assert_eq!(ins.sequence_range(0, 11), "APRTSETMNVEPRF");
    assert!(ins.has_indels() && ins.has_inserts && !ins.has_deletes);

    let del = SeqLoci::new(&mut al, p, seq::parse_reference_sequence(".PRTS.NVEPR."));
    assert_eq!(del.len(), 12);
    assert_eq!(del.at(0), tok::DEL);
    assert_eq!(del.at(5), tok::DEL);
    assert_eq!(del.at(11), tok::DEL);
    assert!(del.has_indels() && !del.has_inserts && del.has_deletes);

    let mut w = SeqLoci::new(&mut al, a, seq::create_aligned("APRT***MNV**EPRF", r));
    assert!(!w.has_inserts && !w.has_deletes && w.has_wildcards());
    w.set_exon_boundary_wildcards(&[2]);
    assert!(!w.has_exon_boundary_wildcards);
    w.set_exon_boundary_wildcards(&[5]);
    assert!(w.has_exon_boundary_wildcards);
}

#[test]
fn hla_sequence_matching() {
    let mut al = Alleles::new();
    let a = al.request("A*01:01:01:01").unwrap();
    let reference = "ABCDEFGHIJKLMNOP";
    let allele_seq = "ABXYEFWZIJKLABOP";
    let s = SeqLoci::new(&mut al, a, seq::create_aligned(allele_seq, reference));
    let loci: Vec<i32> = (0..reference.len()).filter(|i| allele_seq.as_bytes()[*i] != reference.as_bytes()[*i]).map(|i| i as i32).collect();
    let target = |src: &str, loci: &[i32]| -> Vec<Tok> { loci.iter().map(|l| tok::of_byte(src.as_bytes()[*l as usize])).collect() };
    assert_eq!(s.match_type_toks(&target(reference, &loci), &loci), MatchType::Mismatch);
    assert_eq!(s.match_type_toks(&target(allele_seq, &loci), &loci), MatchType::Full);
    let w = SeqLoci::new(&mut al, a, seq::create_aligned("ABXYEFW*********", reference));
    let frag = "ABXYEFWZIJKLABOP";
    assert_eq!(w.match_type_toks(&target(frag, &[5, 6, 7, 8]), &[5, 6, 7, 8]), MatchType::Wild);
    assert_eq!(w.match_type_toks(&target(frag, &[8, 9, 10, 11]), &[8, 9, 10, 11]), MatchType::Wild);
}

// ---------------------------------------------------------------- fragments

#[test]
fn amino_acid_qual_filter() {
    let p = params();
    let pooled = counts(&[(0, "M", 100), (1, "M", 100)]);
    let mut f = empty_fragment("01");
    f.aas = fragment::LociMap::from_pairs(vec![(0, t("M")), (1, t("V"))]);
    let local = counts(&[(0, "M", 10), (1, "M", 10)]);
    fragment::apply_aa_qual_filter(&mut f, &pooled, Some(&local), &p);
    assert_eq!(f.aas.keys().collect::<Vec<_>>(), vec![0]);

    let mut f = empty_fragment("01");
    f.aas = fragment::LociMap::from_pairs(vec![(0, t("M")), (1, t("V"))]);
    let local = counts(&[(0, "M", 10), (1, "M", 9)]);
    fragment::apply_aa_qual_filter(&mut f, &pooled, Some(&local), &p);
    assert_eq!(f.aas.keys().collect::<Vec<_>>(), vec![0, 1]);
}

#[test]
fn nucleotides_from_amino_acids() {
    let s = |x: [Tok; 3]| x.iter().map(|t| tok::to_string(*t)).collect::<Vec<_>>();
    assert_eq!(s(fragment::nucleotides_from_amino_acid(t("X")).unwrap()), vec!["T", "A", "A"]);
    assert_eq!(s(fragment::nucleotides_from_amino_acid(t(".")).unwrap()), vec![".", ".", "."]);
    assert_eq!(s(fragment::nucleotides_from_amino_acid(t("SX")).unwrap()), vec!["A", "G", "CTAA"]);
}

#[test]
fn gene_enrichment() {
    let gc = class1_v37();
    let e = GeneEnrichment::create(&gc.aa_boundaries).unwrap();
    let check = |expected: &[HlaGene], gene: HlaGene, aa: &[i32]| {
        let loci = fragment::expand_indices(aa);
        let bases: Vec<&str> = loci.iter().map(|_| "G").collect();
        let mut f = frag_from_quals("01", gene, vec![gene], &loci, &vec![0; loci.len()], &bases);
        fragment::check_additional_genes(&e, &gc, &mut f);
        let mut got = f.genes.clone();
        got.sort();
        let mut exp = expected.to_vec();
        exp.sort();
        assert_eq!(got, exp, "gene {gene:?} loci {aa:?}");
    };
    check(&[HlaGene::A, HlaGene::B], HlaGene::A, &[337, 348]);
    check(&[HlaGene::A, HlaGene::B], HlaGene::B, &[337, 348]);
    check(&[HlaGene::C], HlaGene::C, &[337, 348]);
    check(&[HlaGene::A], HlaGene::A, &[337, 348, 365]);
    check(&[HlaGene::B], HlaGene::B, &[337, 348, 365]);
}

#[test]
fn nucleotide_quality_filter() {
    let p = params();
    let pooled = counts(&[(0, "A", 100), (1, "A", 100)]);
    let f = frag_from_quals("01", HlaGene::A, vec![HlaGene::A], &[0, 1], &[40, 40], &["A", "T"]);
    let local = counts(&[(0, "A", 10), (1, "A", 10)]);
    let r = fragment::apply_quality_filter(&f, &pooled, &pooled, Some(&local), &p).expect("a new fragment");
    assert_eq!(r.nucs.keys().collect::<Vec<_>>(), vec![0]);
    // both loci below the depth filter (the Java test sets the first locus to 9 and leaves the second empty)
    let local = counts(&[(0, "A", 9)]);
    assert!(fragment::apply_quality_filter(&f, &pooled, &pooled, Some(&local), &p).is_none());
}

#[test]
fn amino_acid_indices() {
    let range = |v: Vec<i32>| (v[0], *v.last().unwrap());
    assert_eq!(range(fragment::calc_amino_acid_indices(0, 0)), (0, -1));
    assert_eq!(range(fragment::calc_amino_acid_indices(0, 1)), (0, -1));
    assert_eq!(range(fragment::calc_amino_acid_indices(0, 2)), (0, 0));
    assert_eq!(range(fragment::calc_amino_acid_indices(1, 2)), (1, 0));
    assert_eq!(range(fragment::calc_amino_acid_indices(2, 2)), (1, 0));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 5)), (1, 1));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 6)), (1, 1));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 7)), (1, 1));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 8)), (1, 2));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 9)), (1, 2));
    assert_eq!(range(fragment::calc_amino_acid_indices(3, 10)), (1, 2));
}

fn create_read(region: (i32, i32), r: &Rec, sc: bool, indels: bool) -> Read {
    let p = params();
    let ctx = ReadCtx { params: &p, ref_bases: None };
    Read::create(&ctx, region, &Arc::new(r.clone()), sc, indels).unwrap()
}

#[test]
fn read_positions() {
    let bases = random_bases(100);
    let r = sam_record("READ_ID_001", "1", 120, &bases[..50], "10S30M10S", "1", 300, false);
    let x = create_read((100, 200), &r, false, false);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end, x.read_index_start, x.read_index_end), (120, 149, 0, 0, 10, 39));
    let x = create_read((100, 200), &r, true, false);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end, x.read_index_start, x.read_index_end), (110, 159, 10, 10, 0, 49));
    let x = create_read((115, 155), &r, true, false);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end, x.read_index_start, x.read_index_end), (115, 155, 5, 6, 5, 45));
    let x = create_read((125, 135), &r, true, false);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end, x.read_index_start, x.read_index_end), (125, 135, 0, 0, 15, 25));
}

#[test]
fn read_low_qual_trimming() {
    let bases = random_bases(100);
    let mut r = sam_record("READ_ID_001", "1", 100, &bases[..50], "50M", "1", 300, false);
    let set = |r: &mut Rec, s: usize, e: usize, q: u8| {
        for i in s..=e {
            r.quals[i] = q;
        }
    };
    set(&mut r, 40, 49, 20);
    let x = create_read((100, 1000), &r, true, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end), (10, 100, 139, 0, 39));
    set(&mut r, 28, 31, 20);
    let x = create_read((100, 1000), &r, false, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end), (10, 100, 139, 0, 39));
    let x = create_read((100, 130), &r, false, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end), (10, 100, 130, 0, 30));
    set(&mut r, 25, 49, 20);
    let x = create_read((100, 130), &r, false, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end), (25, 100, 124, 0, 24));
    r.flags |= 0x10;
    set(&mut r, 0, 49, 37);
    set(&mut r, 0, 9, 20);
    let x = create_read((100, 200), &r, false, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end), (10, 110, 149, 10, 49));

    let mut r = sam_record("READ_ID_001", "1", 120, &bases[..50], "20S30M", "1", 300, false);
    r.flags |= 0x10;
    set(&mut r, 0, 9, 20);
    let x = create_read((100, 200), &r, true, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end, x.soft_clipped_start, x.soft_clipped_end), (10, 110, 149, 10, 49, 10, 0));
    set(&mut r, 0, 29, 20);
    let x = create_read((100, 200), &r, true, false);
    assert_eq!((x.trimmed_bases, x.position_start, x.position_end, x.read_index_start, x.read_index_end, x.soft_clipped_start, x.soft_clipped_end), (30, 130, 149, 30, 49, 0, 0));
}

#[test]
fn adapter_soft_clipped_fragments() {
    let bases = random_bases(100);
    let mut r = sam_record("READ_ID_001", "1", 200, &bases[..80], "10S60M10S", "1", 300, false);
    r.isize = 200;
    let x = create_read((100, 1000), &r, true, true);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end), (190, 269, 10, 10));
    r.isize = 70;
    r.flags |= 0x10;
    let x = create_read((100, 1000), &r, true, true);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end), (200, 269, 0, 10));
    r.flags &= !0x10;
    let x = create_read((100, 1000), &r, true, true);
    assert_eq!((x.position_start, x.position_end, x.soft_clipped_start, x.soft_clipped_end), (190, 259, 10, 0));
}

#[test]
fn fragment_add_nucleotide() {
    let mut f = frag_from_quals("01", HlaGene::A, vec![HlaGene::A], &[1, 2, 3, 6, 7, 8], &[37, 25, 37, 37, 37, 25], &["A", "G", "T", "C", "A", "G"]);
    f.remove_low_qual_bases();
    assert!(!f.contains_nucleotide(2));
    let nuc = |b: &str, q: i8| Nuc { bases: t(b), low_qual: fragment::is_low_qual_byte(q), qual: q };
    f.add_nucleotide(5, nuc("G", 30));
    assert!(f.contains_nucleotide(5));
    f.add_nucleotide(9, nuc("T", 30));
    assert!(f.contains_nucleotide(9));
    assert!(f.validate());
    f.add_nucleotide(0, nuc("A", 37));
    assert!(f.contains_nucleotide(0));
    assert!(f.validate());
}

#[test]
fn fragment_merge() {
    let mut f1 = frag_from_quals("01", HlaGene::A, vec![HlaGene::A], &[1], &[30], &["A"]);
    let f2 = frag_from_quals("01", HlaGene::B, vec![HlaGene::B], &[1], &[30], &["A"]);
    fragment::merge_fragments(&mut f1, &f2);
    assert!(f1.validate());
    assert_eq!(f1.genes.len(), 2);
    assert_eq!(f1.nucs.len(), 1);
    assert_eq!(f1.min_nucleotide_locus(), 1);
    let f2 = frag_from_quals("01", HlaGene::A, vec![HlaGene::A], &[0, 1, 2, 3], &[30; 4], &["A"; 4]);
    fragment::merge_fragments(&mut f1, &f2);
    assert_eq!(f1.genes.len(), 2);
    assert_eq!(f1.nucs.keys().collect::<Vec<_>>(), vec![0, 1, 2, 3]);
    let f2 = frag_from_quals("01", HlaGene::C, vec![HlaGene::C], &[3, 4, 5], &[30; 3], &["A"; 3]);
    fragment::merge_fragments(&mut f1, &f2);
    assert!(f1.validate());
    assert_eq!(f1.genes.len(), 3);
    assert_eq!(f1.nucs.keys().collect::<Vec<_>>(), vec![0, 1, 2, 3, 4, 5]);
}

// ---------------------------------------------------------------- coverage

#[test]
fn allele_coverage_expansion() {
    let gc = class1_v37();
    let mut al = Alleles::new();
    let names = ["A*01:01", "A*02:01", "B*01:01", "B*02:01", "C*01:01", "C*02:01"];
    let ids: Vec<AId> = names.iter().map(|n| al.request_four_digit(n).unwrap()).collect();
    let cov = |i: usize| AlleleCoverage::new(ids[i], 10, 20.0, 0.0);
    let c = ComplexCoverage::create(&gc, &al, (0..6).map(cov).collect());
    assert_eq!(c.unique, 60);
    assert_eq!(c.total, 180);
    let mut c = ComplexCoverage::create(&gc, &al, vec![cov(0), cov(2), cov(4)]);
    c.expand_to_six(&gc, &al);
    assert_eq!(c.alleles.len(), 6);
    assert_eq!(c.unique, 30);
    assert_eq!(c.total, 90);
    for a in &c.alleles {
        assert!((a.total - 15.0).abs() < 0.01);
    }
    let mut c = ComplexCoverage::create(&gc, &al, vec![cov(0), cov(1), cov(4)]);
    c.populate_missing(&gc, &ids);
    assert_eq!(c.alleles.len(), 6);
    assert_eq!(c.total, 90);
    let totals: Vec<f64> = c.alleles.iter().map(|x| x.total).collect();
    assert_eq!(totals, vec![30.0, 30.0, 0.0, 0.0, 30.0, 0.0]);
}

#[test]
fn wildcard_filtering() {
    let mut al = Alleles::new();
    let a1 = al.request_four_digit("A*01:01").unwrap();
    let a2 = al.request_four_digit("A*01:02").unwrap();
    let w3 = al.request_four_digit("A*01:56").unwrap();
    let w4 = al.request_four_digit("A*01:77").unwrap();
    let w5 = al.request_four_digit("A*01:102").unwrap();
    for w in [w3, w4, w5] {
        al.set_has_wildcards(w);
    }
    let mut frags: Vec<Fragment> = (1..=8).map(|i| empty_fragment(&format!("0{i}"))).collect();
    let fa = |i: usize, full: &[AId], wild: &[AId]| FragAlleles::new(&al, i, full, wild);
    let mut fas = vec![
        fa(0, &[w3], &[]),
        fa(1, &[w3, w4, w5], &[]),
        fa(2, &[a1], &[w5]),
        fa(3, &[a2], &[w5]),
        fa(4, &[], &[w3, w4, w5]),
        fa(5, &[], &[w4]),
        fa(6, &[], &[w5]),
        fa(7, &[], &[w4, w5]),
    ];
    let wc = coverage::find_wildcard_alleles(&al, &fas);
    assert_eq!(wc.len(), 3);
    let unsupported = coverage::find_unsupported_wildcards(&al, &fas, &wc);
    assert_eq!(unsupported.len(), 2);
    assert!(unsupported.contains(&w4) && unsupported.contains(&w5));
    coverage::filter_unsupported_wildcard_fragments(&mut fas, &mut frags, &unsupported);
    assert_eq!(fas.len(), 5);
    for id in ["06", "07", "08"] {
        assert!(!fas.iter().any(|x| frags[x.fi].id() == id));
    }
    assert_eq!(frags[5].scope, Scope::WildOnly);
}

#[test]
fn complex_building_two_confirmed_proteins() {
    let mut al = Alleles::new();
    // alleles without a cached group: each is its own group, as in the Java test
    let mk = |al: &mut Alleles, g: &str, p: &str| -> AId {
        let a = al.request_four_digit(&format!("A*{g}:{p}")).unwrap();
        al.v[a as usize].group_allele = None;
        a
    };
    let a11 = mk(&mut al, "1", "1");
    let a12 = mk(&mut al, "1", "2");
    let a13 = mk(&mut al, "1", "3");
    let a21 = mk(&mut al, "2", "1");
    let a22 = mk(&mut al, "2", "2");
    let a23 = mk(&mut al, "2", "3");
    let a31 = mk(&mut al, "3", "1");
    let all = vec![a11, a12, a13, a21, a22, a23, a31];
    let groups = vec![al.group(a11), al.group(a12)];
    let r = ComplexBuilder::build_complexes_by_gene(&al, HlaGene::A, &groups, &all);
    assert_eq!(r[0], vec![a11, a12]);
}

fn test_refdata() -> (LilacConfig, RefData) {
    let config = LilacConfig::for_tests("SAMPLE_ID");
    let rd = RefData::new(&config, GeneSelector::MhcClass1).unwrap();
    (config, rd)
}

#[test]
fn fragment_allele_assignment() {
    let (_, mut rd) = test_refdata();
    let p = params();
    let s1 = seq_of(&mut rd.alleles, "A*01:01", &["A", "B", "C", "D"]);
    let a1 = s1.allele;
    let mut het = std::collections::BTreeMap::new();
    for (i, x) in ["A", "B", "C", "D"].iter().enumerate() {
        het.insert(i as i32, vec![t(x)]);
    }
    let mut mapper = Mapper { aa_het_loci: vec![(HlaGene::A, het)], nuc_het_loci: Vec::new(), ref_nucleotides: Vec::new(), stop_loss: Vec::new(), threads: 1 };
    let mut f1 = Fragment::new(read_record("01"), HlaGene::A, vec![HlaGene::A], Vec::<(i32, Nuc)>::new());
    f1.aas = fragment::LociMap::from_pairs(vec![(0, t("A")), (2, t("C")), (3, t("D"))]);
    let mut frags = vec![f1];
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s1], &[]);
    assert_eq!(fas.len(), 1);
    assert!(fas[0].has_full(a1));

    let mut bhet = std::collections::BTreeMap::new();
    bhet.insert(1, vec![t("B")]);
    bhet.insert(3, vec![t("D")]);
    mapper.aa_het_loci.push((HlaGene::B, bhet));
    let s2 = seq_of(&mut rd.alleles, "B*01:01", &["A", "B", "C", "D"]);
    let mut f2 = Fragment::new(read_record("02"), HlaGene::B, vec![HlaGene::B], Vec::<(i32, Nuc)>::new());
    f2.aas = fragment::LociMap::from_pairs(vec![(0, t("L")), (1, t("B")), (2, t("L")), (3, t("D"))]);
    let mut frags = vec![f2];
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s2], &[]);
    assert_eq!(fas.len(), 1);
    assert!(fas[0].has_full(s2.allele));

    let s3 = seq_of(&mut rd.alleles, "A*01:22", &["*", "B", "C", "*"]);
    let mut f3 = Fragment::new(read_record("03"), HlaGene::A, vec![HlaGene::A], Vec::<(i32, Nuc)>::new());
    f3.aas = fragment::LociMap::from_pairs(vec![(0, t("A")), (1, t("B")), (2, t("C")), (3, t("D"))]);
    let mut frags = vec![f3];
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s3], &[]);
    assert!(fas.is_empty());
    assert_eq!(frags[0].scope, Scope::WildOnly);
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s1, &s3], &[]);
    assert_eq!(fas.len(), 1);
    assert!(fas[0].has_full(a1));
    assert!(fas[0].has_wild(s3.allele));
}

#[test]
fn fragment_allele_low_base_quals() {
    let (_, mut rd) = test_refdata();
    let p = params();
    let s1 = seq_of(&mut rd.alleles, "A*01:01", &["A", "P", "C", "D"]);
    let mut het = std::collections::BTreeMap::new();
    for (i, x) in ["A", "P", "C", "D"].iter().enumerate() {
        het.insert(i as i32, vec![t(x)]);
    }
    let mapper = Mapper { aa_het_loci: vec![(HlaGene::A, het)], nuc_het_loci: Vec::new(), ref_nucleotides: Vec::new(), stop_loss: Vec::new(), threads: 1 };
    let sequence = "GCTCCCTGCGAC";
    let loci: Vec<i32> = (0..sequence.len() as i32).collect();
    let bases: Vec<String> = sequence.chars().map(|c| c.to_string()).collect();
    let mut quals: Vec<i8> = vec![30; sequence.len()];
    quals[1] = 2;
    quals[4] = 2;
    quals[11] = 2;
    let mut f = frag_from_quals("01", HlaGene::A, vec![HlaGene::A], &loci, &quals, &bases.iter().map(|x| x.as_str()).collect::<Vec<_>>());
    assert!(f.validate());
    f.remove_low_qual_bases();
    f.build_amino_acids();
    assert_eq!(f.nucs.len(), 9);
    assert_eq!(f.aas.len(), 1);
    let mut frags = vec![f];
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s1], &[]);
    assert_eq!(fas.len(), 1);
    assert!(fas[0].has_full(s1.allele));
    // the second raw base changed so its codon no longer gives the amino acid: still a match
    let second = frags[0].raw.keys().nth(1).unwrap();
    let mut n = frags[0].raw.get(second).unwrap();
    n.bases = t("A");
    frags[0].raw.insert(second, n);
    let fas = mapper.create_fragment_alleles(&rd, &p, &mut frags, &[0], &[&s1], &[]);
    assert_eq!(fas.len(), 1);
    assert!(fas[0].has_full(s1.allele));
}

#[test]
fn fragment_allele_matrix_coverage() {
    let mut al = Alleles::new();
    let a1 = al.request_four_digit("A*01:01").unwrap();
    let a2 = al.request_four_digit("B*01:01").unwrap();
    let a3 = al.request_four_digit("C*01:01").unwrap();
    let fas = vec![
        FragAlleles::new(&al, 0, &[a1], &[]),
        FragAlleles::new(&al, 1, &[a2], &[]),
        FragAlleles::new(&al, 2, &[a3], &[]),
        FragAlleles::new(&al, 3, &[a2], &[a1, a3]),
        FragAlleles::new(&al, 4, &[], &[a1, a2, a3]),
    ];
    let m = crate::complexes::Matrix::new(&fas, &[a1, a2, a3]);
    let c = m.create(&[a1, a2, a3]);
    assert_eq!(c.len(), 3);
    for x in &c {
        assert!((x.total - 1.67).abs() < 0.01);
        assert_eq!(x.unique, 1);
    }
    assert!((c[0].wild - 0.67).abs() < 0.01);
    assert!((c[1].wild - 0.33).abs() < 0.01);
    assert!((c[2].wild - 0.67).abs() < 0.01);
}

#[test]
fn homozygous_counts() {
    let gc = class1_v37();
    let mut al = Alleles::new();
    let mut cov = |names: &[&str]| -> ComplexCoverage {
        let l: Vec<AlleleCoverage> = names.iter().map(|n| AlleleCoverage::new(al.request_four_digit(n).unwrap(), 10, 0.0, 0.0)).collect();
        ComplexCoverage::create(&gc, &al, l)
    };
    assert_eq!(cov(&["A*01:02", "B*01:01", "C*03:01"]).homozygous_count, 3);
    assert_eq!(cov(&["A*01:02", "A*01:03", "B*01:01", "B*01:02", "C*03:01"]).homozygous_count, 1);
    assert_eq!(cov(&["A*01:02", "A*01:03", "B*01:01", "B*01:02", "C*03:01", "C*03:02"]).homozygous_count, 0);
}

#[test]
fn coverage_rank_sorting() {
    let mut al = Alleles::new();
    let complexes: Vec<(Vec<AId>, f64)> = vec![
        (["A*01:201", "B*01:01", "C*01:01"].iter().map(|n| al.request_four_digit(n).unwrap()).collect(), 34.0),
        (["A*02:01", "B*02:01", "C*02:01"].iter().map(|n| al.request_four_digit(n).unwrap()).collect(), 33.0),
        (["A*03:01", "B*03:01", "C*03:01"].iter().map(|n| al.request_four_digit(n).unwrap()).collect(), 35.0),
        (["A*01:78", "B*01:02", "C*01:02"].iter().map(|n| al.request_four_digit(n).unwrap()).collect(), 34.0),
    ];
    let (config, mut rd) = test_refdata();
    rd.alleles = al;
    let ranker = crate::complexes::Ranker::new(&rd, 0.0, 3, config.params.freq_score_penalty);
    let items: Vec<crate::complexes::Scored> = complexes
        .iter()
        .enumerate()
        .map(|(i, c)| crate::complexes::Scored {
            index: i as u32,
            totals: Default::default(),
            score: c.1,
            cohort_total: 0.0,
            cohort_penalty: 0.0,
            complexity: 0.0,
            complexity_penalty: 0.0,
            recovered: 0,
            wildcards: 0,
        })
        .collect();
    let order = ranker.rank(&items, &|s: &crate::complexes::Scored| complexes[s.index as usize].0.clone());
    assert_eq!(order, vec![2, 3, 0, 1]);
}

fn complexity(seqs: &[(&str, &[&str])]) -> i32 {
    let mut al = Alleles::new();
    let mut exons = std::collections::HashMap::new();
    let mut ids = Vec::new();
    for (n, s) in seqs {
        let sl = seq_of(&mut al, n, s);
        exons.insert(sl.allele, seq::exon_sequences(&[1, 2, 100], &sl));
        ids.push(sl.allele);
    }
    ids.sort_by(|a, b| al.name(*a).cmp(al.name(*b)));
    let mut cx = crate::complexes::Complexity::new();
    crate::complexes::complexity_count(&mut cx, &exons, &al, &ids)
}

#[test]
fn solution_complexity() {
    let aaa: &[&str] = &["A", "A", "A", "A", "A", "A"];
    assert_eq!(complexity(&[("A*01:02", aaa)]), 3);
    assert_eq!(complexity(&[("A*01:02", aaa), ("A*01:03", aaa), ("A*02:02", aaa)]), 3);
    assert_eq!(complexity(&[("A*01:02", aaa), ("A*01:03", &["A", "B", "A", "A", "A", "A"]), ("A*02:02", &["A", "A", "B", "A", "A", "A"])]), 5);
    assert_eq!(complexity(&[("A*01:02", aaa), ("A*01:03", &["A", "B", "A", "A", "B", "A"]), ("A*02:02", &["A", "A", "B", "A", "*", "A"])]), 5);
    assert_eq!(complexity(&[("A*01:02", aaa), ("A*01:03", &["A", "B", "A", "A", "B", "A"]), ("A*02:02", &["A", "A", "B", "A", "*", "A", "C"])]), 7);
}

// ---------------------------------------------------------------- evidence

#[test]
fn filter_sequences_by_min_support() {
    let p = params();
    let mut al = Alleles::new();
    let s1 = seq_of(&mut al, "A*01:01", &["M", "M"]);
    let s2 = seq_of(&mut al, "A*01:01", &["M", "V"]);
    let pooled = counts(&[(0, "M", 100), (1, "M", 100)]);
    let local = counts(&[(0, "M", 10), (1, "M", 10)]);
    let c = evidence::filter_by_min_support(&p, vec![&s1, &s2], &pooled, &[], Some(&local));
    assert_eq!(c.len(), 1);
    assert_eq!(c[0].sequence_str(), "MM");
    let local = counts(&[(0, "M", 10), (1, "M", 9)]);
    assert_eq!(evidence::filter_by_min_support(&p, vec![&s1, &s2], &pooled, &[], Some(&local)).len(), 2);
}

#[test]
fn filter_candidates_by_phased_evidence() {
    let p = params();
    let mut al = Alleles::new();
    let s1 = seq_of(&mut al, "A*01:01", &["M", "L", "V"]);
    let s2 = seq_of(&mut al, "A*01:01", &["M", "L", "L"]);
    let mut ev = vec![PhasedEvidence::new(vec![0, 2], vec![("MV".to_string(), 0)])];
    let local = counts(&[(0, "M", 10), (1, "M", 10), (2, "M", 10)]);
    let c = evidence::filter_candidates(&p, vec![&s1, &s2], &mut ev, Some(&local));
    assert_eq!(c.len(), 1);
    assert_eq!(c[0].sequence_str(), "MLV");
    let local = counts(&[(0, "M", 10), (2, "M", 9)]);
    let mut ev = vec![PhasedEvidence::new(vec![0, 2], vec![("MV".to_string(), 0)])];
    assert_eq!(evidence::filter_candidates(&p, vec![&s1, &s2], &mut ev, Some(&local)).len(), 2);
}

fn ev(loci: &[i32], e: &[(&str, i32)]) -> PhasedEvidence {
    PhasedEvidence::new(loci.to_vec(), e.iter().map(|x| (x.0.to_string(), x.1)).collect())
}

#[test]
fn combine_evidence() {
    let left = ev(&[344], &[("T", 13), ("S", 71)]);
    let right = ev(&[348], &[("S", 20), ("C", 20)]);
    let combined = ev(&[344, 348], &[("SS", 9), ("SC", 4)]);
    assert!(!evidence::can_combine3(&left, &combined, &right));

    let left = ev(&[1, 3], &[("AM", 36), ("LM", 16), ("RT", 11), ("RM", 35)]);
    let right = ev(&[1, 3, 7], &[("AMT", 31), ("AMA", 1), ("LMT", 15), ("RTT", 9), ("RTA", 1), ("RMA", 16)]);
    assert!(evidence::can_combine(&left, &right));

    let left = ev(&[3, 7], &[("AR", 2), ("TT", 3)]);
    assert!(evidence::can_combine(&left, &ev(&[3, 7, 10], &[("ARS", 2), ("TTH", 3)])));
    assert!(!evidence::can_combine(&left, &ev(&[3, 7, 10], &[("ARS", 2)])));
    assert!(evidence::can_combine(&ev(&[2, 3, 7], &[("HAR", 2), ("ATT", 3)]), &left));
    let right = ev(&[7, 10], &[("RS", 2), ("TS", 3)]);
    assert!(evidence::can_combine(&left, &right));
    assert!(!evidence::can_combine(&right, &left));
    assert!(!evidence::can_combine(&left, &ev(&[7, 1], &[("SS", 2), ("TS", 3)])));
    assert!(!evidence::can_combine(&ev(&[3, 7], &[("AR", 2), ("TR", 3)]), &ev(&[9, 10], &[("TS", 2), ("QR", 3)])));
}

#[test]
fn inconsistent_evidence() {
    let mut al = Alleles::new();
    let victim = ev(&[0, 1, 3], &[("CAT", 4), ("ATC", 5)]);
    let cat = seq_loci(&mut al, "A*01:01", "CART");
    let atc = seq_loci(&mut al, "A*01:02", "ATRC");
    let wild_atc = seq_loci(&mut al, "A*01:03", "*TRC");
    let wild = seq_loci(&mut al, "A*01:04", "****");
    assert_eq!(victim.inconsistent_evidence(&[&cat, &atc]).total(), 0);
    let m = victim.inconsistent_evidence(&[&atc]);
    assert_eq!(m.evidence.len(), 1);
    assert_eq!(m.evidence[0].0, "CAT");
    let m = victim.inconsistent_evidence(&[&wild_atc]);
    assert_eq!(m.evidence.len(), 1);
    assert_eq!(m.evidence[0].0, "CAT");
    assert_eq!(victim.inconsistent_evidence(&[&wild]).total(), 0);
}

// ---------------------------------------------------------------- misc

#[test]
fn nucleotide_locus() {
    for v38 in [false, true] {
        let tr = refdata::load_hla_transcripts(v38, Some(GeneSelector::MhcClass1));
        let get = |g: HlaGene| tr.iter().find(|x| x.0 == g).unwrap().1.clone();
        let all = || tr.iter().map(|x| &x.1);
        let (a, b, c) = (get(HlaGene::A), get(HlaGene::B), get(HlaGene::C));
        assert_eq!(refdata::nucleotide_locus(all(), a.coding_start), 0);
        assert_eq!(refdata::nucleotide_locus(all(), a.coding_end), 1097);
        assert_eq!(refdata::nucleotide_locus(all(), b.coding_end), 0);
        assert_eq!(refdata::nucleotide_locus(all(), b.coding_start), 1088);
        assert_eq!(refdata::nucleotide_locus(all(), c.coding_end), 0);
        assert_eq!(refdata::nucleotide_locus(all(), c.coding_start), 1100);
    }
}

#[test]
fn unmatched_haplotype() {
    let mut al = Alleles::new();
    let counts = SequenceCount::from_counts(&[vec![("C", 4), ("A", 4)], vec![("A", 4), ("T", 4)], vec![("R", 8)], vec![("T", 4), ("C", 4)]]);
    let victim = ev(&[0, 1, 3], &[("CAT", 4), ("ATC", 5)]);
    let cat = seq_loci(&mut al, "A*01:01", "CART");
    let atc = seq_loci(&mut al, "A*01:02", "ATRC");
    let wild_atc = seq_loci(&mut al, "A*01:03", "*TRC");
    let wild = seq_loci(&mut al, "A*01:04", "****");
    assert!(qc::unmatched_haplotype(&victim, 0, &[&cat, &atc], &counts, Some(&[])).is_empty());
    let m = qc::unmatched_haplotype(&victim, 0, &[&atc], &counts, Some(&[]));
    assert_eq!(m.len(), 1);
    assert_eq!(m[0].haplotype, b"CART");
    assert!(qc::unmatched_haplotype(&victim, 5, &[&atc], &counts, Some(&[])).is_empty());
    let m = qc::unmatched_haplotype(&victim, 0, &[&wild_atc], &counts, Some(&[]));
    assert_eq!(m.len(), 1);
    assert_eq!(m[0].haplotype, b"CART");
    assert!(qc::unmatched_haplotype(&victim, 0, &[&wild], &counts, Some(&[])).is_empty());
}

#[test]
fn unmatched_amino_acid() {
    let mut al = Alleles::new();
    let w1 = seq_loci(&mut al, "A*01:01", "ABCDEFGHIJKLMNO");
    let w2 = seq_loci(&mut al, "A*01:02", "ABCDEFGHIJKLM**");
    let w3 = seq_loci(&mut al, "A*01:02", "A**DEFGHIJKLMNN");
    let warn = (0.01f64 * 1000.0) as u32;
    let mut maps: Vec<Vec<(String, u32)>> = (0..w1.len()).map(|i| vec![(tok::to_string(w1.at(i)), 10)]).collect();
    maps[1].push(("C".into(), warn));
    maps[14].push(("C".into(), warn));
    maps[7].push(("C".into(), warn));
    maps[3].push(("C".into(), warn - 1));
    maps[0].push(("B".into(), warn));
    maps[12].push(("B".into(), warn));
    let refs: Vec<Vec<(&str, u32)>> = maps.iter().map(|m| m.iter().map(|x| (x.0.as_str(), x.1)).collect()).collect();
    let counts = SequenceCount::from_counts(&refs);
    let unmatched = vec![Haplotype { start: 5, end: 10, phased_fragments: 10, haplotype: b"ABCDEF".to_vec(), matching_fragments: 0 }];
    let q = AminoAcidQC::create(&[&w1, &w2, &w3], Some(&[]), &counts, &unmatched);
    assert_eq!(q.unused, 3);
    assert_eq!(q.max_frags, 10);
}

#[test]
fn somatic_coding_count_repeat_allele() {
    let mut al = Alleles::new();
    let names = ["A*01:01", "A*01:01", "B*01:01", "B*01:02", "C*01:01", "C*01:03"];
    let ids: Vec<AId> = names.iter().map(|n| al.request_four_digit(n).unwrap()).collect();
    let mut c: Vec<SomaticCodingCount> = ids.iter().map(|a| SomaticCodingCount::new(*a)).collect();
    assert_eq!(c.len(), 6);
    SomaticCodingCount::add_variant(&mut c, &al, true, CodingEffect::Missense, &[ids[0], ids[2]]);
    SomaticCodingCount::add_variant(&mut c, &al, false, CodingEffect::Missense, &[ids[2]]);
    assert!((c[0].total() - 0.5).abs() < 0.001);
    assert!((c[0].inframe_indel - 0.5).abs() < 0.001);
    assert!((c[1].total()).abs() < 0.001);
    assert!((c[2].total() - 1.5).abs() < 0.001);
    assert!((c[2].inframe_indel - 0.5).abs() < 0.001);
    assert!((c[2].missense - 1.0).abs() < 0.001);
}

// ---------------------------------------------------------------- application

struct AppTest {
    config: LilacConfig,
    rd: RefData,
    counter: usize,
}

impl AppTest {
    /// createLilacApp and loadTestReferenceData: the test sequences, frequencies of 0.1 for the common alleles
    fn new() -> AppTest {
        let mut config = LilacConfig::for_tests("SAMPLE_ID");
        // disableCoverageThresholds
        config.params.fail_low_coverage_threshold = 500;
        config.params.warn_low_coverage_threshold = 500.0;
        let mut rd = RefData::new(&config, GeneSelector::MhcClass1).unwrap();
        let nuc = hmf_common::java_text::lines_of(include_str!("../testdata/test_allele_nucleotides.csv"));
        let aa = hmf_common::java_text::lines_of(include_str!("../testdata/test_allele_amino_acids.csv"));
        let mut template = None;
        rd.nuc_seqs = rd.load_sequence_lines(&nuc, false, &[], &[], &mut template);
        rd.aa_seqs = rd.load_sequence_lines(&aa, true, &[], &[], &mut template);
        for name in ["A*01:01", "A*02:01", "A*03:01", "B*07:02", "B*08:01", "C*01:02", "C*02:02"] {
            if rd.aa_seqs.iter().any(|s| rd.alleles.name(s.allele) == name) {
                rd.freqs.put(name, 0.1);
            }
        }
        AppTest { config, rd, counter: 0 }
    }

    fn fragments(&mut self, allele: &str, n: usize, start: i32, end: i32, length: i32, gap: i32) -> Vec<Fragment> {
        let mut out = Vec::new();
        let Some(s) = self.rd.nuc_seqs.iter().find(|s| self.rd.alleles.name(s.allele) == allele) else { return out };
        let seq = s.sequence_str();
        let gene = self.rd.alleles.gene(s.allele);
        let mut locus = start;
        for _ in 0..n {
            if locus + length > end {
                locus = start;
            }
            self.counter += 1;
            out.push(test_fragment(&format!("READ_{:03}", self.counter), gene, &seq, locus, locus + length));
            locus += gap;
        }
        out
    }

    fn run(&mut self, frags: Vec<Fragment>, stop_loss: Vec<u64>) -> crate::app::GroupResult {
        crate::app::run_gene_group(&self.config, &mut self.rd, None, crate::app::Fragments::Given { reference: frags, tumor: Vec::new(), stop_loss }, None).unwrap()
    }

    fn names(&self, ids: &[AId]) -> Vec<String> {
        ids.iter().map(|a| self.rd.alleles.name(*a).to_string()).collect()
    }
}

#[test]
fn app_basic() {
    let mut a = AppTest::new();
    let mut frags = Vec::new();
    frags.extend(a.fragments("A*01:01:01", 50, 20, 1050, 150, 20));
    for n in ["A*02:01:01", "B*07:02:01", "B*08:01:01", "C*01:02:01", "C*02:02:01"] {
        frags.extend(a.fragments(n, 200, 20, 1050, 150, 20));
    }
    let r = a.run(frags, Vec::new());
    assert_eq!(r.ranked.len(), 1);
    assert_eq!(r.ranked[0].homozygous_count, 0);
    assert!((r.ranked[0].cohort_frequency_total - -6.0).abs() < 0.01);
    let names = a.names(&r.reference_alleles);
    assert_eq!(names.len(), 6);
    for n in ["A*01:01", "A*02:01", "B*07:02", "B*08:01", "C*01:02", "C*02:02"] {
        assert!(names.contains(&n.to_string()), "{n} in {names:?}");
    }
    let qc = r.qc.unwrap();
    assert!(qc.status.contains(&qc::QcStatus::Pass));
    assert!(qc.hla_y_allele.is_none());
    assert_eq!(qc.aa_qc.unused, 0);
    assert!(qc.haplotype_qc.unmatched.is_empty());
    assert_eq!(qc.coverage_qc.total_fragments, 200 * 6 - 150);
    assert_eq!(qc.coverage_qc.fitted, 200 * 6 - 150);
}

#[test]
fn app_known_stop_loss() {
    let mut a = AppTest::new();
    let mut frags = Vec::new();
    for n in ["A*01:01:01", "A*02:01:01", "B*07:02:01", "B*08:01:01", "C*01:02:01", "C*02:02:01"] {
        frags.extend(a.fragments(n, 50, 20, 1050, 150, 20));
    }
    let stop = a.fragments(refdata::STOP_LOSS_ON_C_ALLELE, 20, 900, 1090, 150, 20);
    let uids: Vec<u64> = stop.iter().map(|f| f.uid).collect();
    frags.extend(stop);
    frags.extend(a.fragments(refdata::STOP_LOSS_ON_C_ALLELE, 50, 20, 1050, 150, 20));
    let r = a.run(frags, uids);
    assert_eq!(r.ranked.len(), 1);
    assert_eq!(r.ranked[0].homozygous_count, 0);
    assert!((r.ranked[0].cohort_frequency_total - -9.0).abs() < 0.01);
    let names = a.names(&r.reference_alleles);
    for n in ["A*01:01", "A*02:01", "B*07:02", "B*08:01", "C*02:02", "C*04:09N"] {
        assert!(names.contains(&n.to_string()), "{n} in {names:?}");
    }
    let qc = r.qc.unwrap();
    assert!(qc.status.contains(&qc::QcStatus::WarnUnmatchedHaplotype));
    assert!(qc.hla_y_allele.is_none());
    assert_eq!(qc.aa_qc.unused, 0);
    assert_eq!(qc.haplotype_qc.unmatched.len(), 5);
    assert_eq!(qc.coverage_qc.total_fragments, 370);
    assert_eq!(qc.coverage_qc.fitted, 338);
}

#[test]
fn app_wildcards() {
    let mut a = AppTest::new();
    let mut frags = Vec::new();
    frags.extend(a.fragments("A*01:01:01", 10, 20, 1050, 150, 20));
    frags.extend(a.fragments("A*02:01:01", 150, 20, 1050, 150, 20));
    frags.extend(a.fragments("A*03:01:01", 300, 20, 1050, 150, 20));
    frags.extend(a.fragments("B*07:02:01", 10, 20, 1050, 150, 20));
    frags.extend(a.fragments("B*08:01:01", 150, 20, 1050, 150, 20));
    frags.extend(a.fragments("B*38:58", 300, 20, 1050, 150, 20));
    frags.extend(a.fragments("C*01:02:01", 150, 20, 1050, 150, 20));
    let r = a.run(frags, Vec::new());
    assert_eq!(r.ranked.len(), 1);
    assert_eq!(r.ranked[0].homozygous_count, 1);
    assert_eq!(r.ranked[0].recovered_count, 0);
    assert!((r.ranked[0].cohort_frequency_total - -9.0).abs() < 0.01);
    let names = a.names(&r.reference_alleles);
    assert_eq!(names.len(), 6);
    for n in ["A*02:01", "A*03:01", "B*08:01", "B*38:58", "C*01:02"] {
        assert!(names.contains(&n.to_string()), "{n} in {names:?}");
    }
    let qc = r.qc.unwrap();
    assert!(qc.status.contains(&qc::QcStatus::Pass));
    assert_eq!(qc.aa_qc.unused, 0);
    assert!(!qc.haplotype_qc.unmatched.is_empty());
    assert!(qc.coverage_qc.percent_wildcard > 0.0);
}
