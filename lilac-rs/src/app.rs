//! LilacApplication: each gene group in turn, from the reads to the typing, tumor coverage and QC.

use std::collections::{BTreeMap, HashSet};

use anyhow::Result;

use crate::bam::{Bam, ReadCtx};
use crate::complexes::{self, apply_score, Product, Ranker, Scored, Totals};
use crate::config::{LilacConfig, WriteType};
use crate::constants::*;
use crate::coverage::{self, calc_protein_coverage, ComplexBuilder, ComplexCoverage, FragAlleles, Mapper};
use crate::evidence::{self, PhasedEvidence};
use crate::fragment::{self, check_additional_genes, AaPipeline, Factory, Fragment, RawCounts, Scope};
use crate::hla::{AId, GeneSelector, HlaGene};
use crate::hlay::{HlaY, Source};
use crate::output::{fragment_line, GroupOut, Writers};
use crate::qc::{AminoAcidQC, BamQC, CoverageQC, HaplotypeQC, LilacQC, SolutionSummary, SomaticCodingCount, SomaticVariantQC};
use crate::reader::GeneReader;
use crate::refdata::{RefBases, RefData, SequenceFiles};
use crate::seq::SeqLoci;
use crate::tok::Tok;
use crate::variant::{self, VariantAnnotation};

/// A source of fragments for one gene group: the BAM files, or fragments given by a test
pub enum Fragments {
    Bams,
    #[allow(dead_code)]
    Given { reference: Vec<Fragment>, tumor: Vec<Fragment>, stop_loss: Vec<u64> },
}

/// What a gene group produced (for the tests)
#[derive(Default)]
pub struct GroupResult {
    pub ranked: Vec<ComplexCoverage>,
    pub summary_lines: Vec<String>,
    pub qc_line: Option<String>,
    /// the reference coverage's alleles of the solution summary
    pub reference_alleles: Vec<AId>,
    pub qc: Option<LilacQC>,
}

pub fn run(config: &LilacConfig) -> Result<()> {
    let started = std::time::Instant::now();
    crate::info!(
        "sample({}) inputs: referenceBam({}) tumorBam({}) somaticVCF({}) geneCopyNumber({}) rnaBam({})",
        config.sample,
        !config.reference_bam.is_empty(),
        !config.tumor_bam.is_empty(),
        !config.somatic_variants_file.is_empty(),
        !config.copy_number_file.is_empty(),
        !config.rna_bam.is_empty()
    );
    if config.run_validation {
        crate::info!("running validation routines");
    }
    if !config.actual_alleles.is_empty() {
        crate::info!("expected alleles: {}", config.actual_alleles.join(", "));
    }
    if !config.restricted_alleles.is_empty() {
        crate::info!("restricted alleles: {}", config.restricted_alleles.join(", "));
    }
    let mut writers = Writers::new(config)?;
    let files = SequenceFiles::load(&config.resource_dir)?;
    let ref_bases = if !config.ref_genome.is_empty() { Some(RefBases::load(&config.ref_genome, config.v38)?) } else { None };
    let run_group = |genes: GeneSelector, actual: &[String], restricted: &[String]| -> Result<(GroupOut, Vec<String>, Vec<String>)> {
        crate::info!("solving for genes({})", genes.name());
        let mut rd = RefData::new(config, genes)?;
        rd.load(&files, config, actual, restricted)?;
        let a: Vec<String> = rd.actual_alleles.iter().map(|x| rd.alleles.name(*x).to_string()).collect();
        let r: Vec<String> = rd.restricted_alleles.iter().map(|x| rd.alleles.name(*x).to_string()).collect();
        let mut out = GroupOut::default();
        run_gene_group(config, &mut rd, ref_bases.as_ref(), Fragments::Bams, Some(&mut out))?;
        Ok((out, a, r))
    };
    // the config's allele lists: as given until the first group's reference data rebuilds them at four digits,
    // so with them the groups run one after the other; otherwise in parallel, written in group order
    let sequential = threads(config) < 2 || config.genes.len() < 2 || !config.actual_alleles.is_empty() || !config.restricted_alleles.is_empty();
    if sequential {
        let mut actual = config.actual_alleles.clone();
        let mut restricted = config.restricted_alleles.clone();
        for genes in &config.genes {
            let (out, a, r) = run_group(*genes, &actual, &restricted)?;
            actual = a;
            restricted = r;
            writers.append(out)?;
        }
    } else {
        let next = std::sync::atomic::AtomicUsize::new(0);
        let n = config.genes.len();
        let mut results: Vec<Option<Result<GroupOut>>> = (0..n).map(|_| None).collect();
        let slots: Vec<std::sync::Mutex<Option<Result<GroupOut>>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|s| {
            for _ in 0..threads(config).min(n) {
                s.spawn(|| loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                    if i >= n {
                        break;
                    }
                    let r = run_group(config.genes[i], &[], &[]).map(|x| x.0);
                    *slots[i].lock().unwrap() = Some(r);
                });
            }
        });
        for (i, s) in slots.into_iter().enumerate() {
            results[i] = s.into_inner().unwrap();
        }
        for r in results {
            writers.append(r.unwrap()?)?;
        }
    }
    writers.close()?;
    crate::info!("Lilac completed, mins({})", hmf_common::fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(())
}

fn threads(config: &LilacConfig) -> usize {
    config.threads.max(1) as usize
}

/// hasSufficientGeneDepth
fn has_sufficient_depth(config: &LilacConfig, depth: &[(HlaGene, Vec<i32>)]) -> bool {
    let mut total = 0;
    let mut per: Vec<String> = Vec::new();
    for (g, d) in depth {
        if *g == HlaGene::Drb3 || *g == HlaGene::Drb4 {
            continue;
        }
        let n = d.iter().filter(|x| (**x as f64) < WARN_LOW_COVERAGE_DEPTH).count() as i32;
        per.push(format!("{}={n}", g.long_name()));
        total += n;
    }
    if total >= config.fatal_total_low_coverage_positions {
        crate::warn!(
            "exiting due to too low coverage: bases with <{} coverage per gene({{{}}}) total({total}) exceeds threshold({})",
            WARN_LOW_COVERAGE_DEPTH,
            per.join(", "),
            config.fatal_total_low_coverage_positions
        );
        return false;
    }
    true
}

/// writeFailedSampleFileOutputs: a FAIL QC line, no alleles
fn write_failed(config: &LilacConfig, rd: &RefData, depth: &[(HlaGene, Vec<i32>)], writers: Option<&mut GroupOut>) {
    let Some(w) = writers else { return };
    if config.output_dir.is_empty() {
        return;
    }
    crate::info!("writing failed-sample output to {}", config.output_dir);
    let gc = &rd.gene_cache;
    let counts: Vec<(HlaGene, i32)> = depth.iter().map(|x| (x.0, 0)).collect();
    let qc = LilacQC::new(
        &config.params,
        rd.genes,
        0.0,
        String::new(),
        None,
        None,
        AminoAcidQC { unused: 0, max_frags: 0 },
        BamQC::new(gc, 0, 0, 0, depth),
        CoverageQC::new(counts, 0, 0, 0, 0, 0, 0, 0),
        HaplotypeQC { unused: 0, max_frags: 0, pon: 0, unmatched: Vec::new() },
        SomaticVariantQC { variant_count: 0, allele_count: 0.0 },
    );
    w.qc.push(qc.body(gc).join("\t"));
}

/// extractHeterozygousLociSequences: per gene the heterozygous amino acids (and those of extra sequences),
/// shared over HLA-A, -B and -C up to locus 337 and over A and B to 348, or over the DRB genes
fn het_aa_loci(rd: &RefData, counts: &[(HlaGene, crate::seq::SequenceCount)], extra: &[&SeqLoci]) -> Vec<(HlaGene, BTreeMap<i32, Vec<Tok>>)> {
    let mut m: Vec<(HlaGene, BTreeMap<i32, Vec<Tok>>)> = counts
        .iter()
        .map(|(g, c)| {
            let ex: Vec<&SeqLoci> = extra.iter().copied().filter(|s| rd.alleles.gene(s.allele) == *g).collect();
            (*g, c.het_loci_sequences(&ex))
        })
        .collect();
    let has = |m: &Vec<(HlaGene, BTreeMap<i32, Vec<Tok>>)>, g: HlaGene| m.iter().any(|x| x.0 == g);
    let union = |a: &[Tok], b: &[Tok]| -> Vec<Tok> {
        let mut v = a.to_vec();
        for t in b {
            if !v.contains(t) {
                v.push(*t);
            }
        }
        v
    };
    if has(&m, HlaGene::A) {
        let get = |m: &Vec<(HlaGene, BTreeMap<i32, Vec<Tok>>)>, g: HlaGene, l: i32| -> Vec<Tok> {
            m.iter().find(|x| x.0 == g).and_then(|x| x.1.get(&l).cloned()).unwrap_or_default()
        };
        for locus in 0..=348 {
            let (a, b, c) = (get(&m, HlaGene::A, locus), get(&m, HlaGene::B, locus), get(&m, HlaGene::C, locus));
            if locus <= 337 {
                if a.is_empty() && b.is_empty() && c.is_empty() {
                    continue;
                }
                let all = union(&union(&a, &b), &c);
                for x in m.iter_mut() {
                    if matches!(x.0, HlaGene::A | HlaGene::B | HlaGene::C) {
                        x.1.insert(locus, all.clone());
                    }
                }
            } else {
                if a.is_empty() && b.is_empty() {
                    continue;
                }
                let ab = union(&a, &b);
                for x in m.iter_mut() {
                    if matches!(x.0, HlaGene::A | HlaGene::B) {
                        x.1.insert(locus, ab.clone());
                    }
                }
            }
        }
    } else if has(&m, HlaGene::Drb3) {
        let loci: std::collections::BTreeSet<i32> = m.iter().flat_map(|x| x.1.keys().copied()).collect();
        for l in loci {
            let mut all: Vec<Tok> = Vec::new();
            for x in &m {
                if let Some(v) = x.1.get(&l) {
                    all = union(&all, v);
                }
            }
            for x in m.iter_mut() {
                x.1.insert(l, all.clone());
            }
        }
    }
    m
}

/// checkDownsampleRefFragmentAlleles: the reference fragment alleles used to rank complexes (indices)
fn downsample(config: &LilacConfig, fas: &[FragAlleles], frags: &[Fragment]) -> Vec<usize> {
    let n = fas.len();
    let max = config.max_ref_fragments;
    if max == 0 || n <= max as usize {
        return (0..n).collect();
    }
    if n <= max as usize * 2 {
        let mut idx: Vec<usize> = (0..n).collect();
        hmf_common::java_sort::java_sort(&mut idx, |a, b| hmf_common::java_text::compare(frags[fas[*a].fi].id(), frags[fas[*b].fi].id()));
        idx.truncate(max as usize);
        return idx;
    }
    let nth = (n as f64 / max as f64).floor() as i64;
    let mut out = Vec::new();
    let mut counter = 0i64;
    for i in 0..n {
        counter += 1;
        if counter == nth {
            out.push(i);
            if out.len() >= max as usize {
                break;
            }
            counter = 0;
        }
    }
    out
}

/// The per-gene work of the gene tasks: phasing fragments, candidates and phased evidence
struct GeneTaskResult {
    candidates: Vec<AId>,
    evidence: Vec<PhasedEvidence>,
}

#[allow(clippy::too_many_arguments)]
fn gene_task(config: &LilacConfig, rd: &RefData, ctx: &crate::refdata::HlaContext, frags: &[Fragment], raw_aa: Option<&crate::seq::SequenceCount>) -> GeneTaskResult {
    let params = &config.params;
    let unphased = evidence::unphased_candidates(params, &rd.alleles, ctx, frags, &rd.aa_seqs, &rd.nuc_seqs, &rd.common_alleles, raw_aa);
    let mut ev = evidence::phased_evidence(config, params, ctx.gene, &ctx.expected, frags);
    let set: HashSet<AId> = unphased.iter().copied().collect();
    crate::debug!("gene({}) determining phased candidate set", ctx.gene.long_name());
    let unphased_seqs: Vec<&SeqLoci> = rd.aa_seqs.iter().filter(|s| set.contains(&rd.alleles.four_digit(s.allele))).collect();
    let phased = evidence::filter_candidates(params, unphased_seqs, &mut ev, raw_aa);
    let alleles: Vec<AId> = phased.iter().map(|s| s.allele).collect();
    crate::debug!("gene({}) has {} candidates after phasing: {}", ctx.gene.long_name(), alleles.len(), rd.alleles.join_max(&alleles, 100));
    GeneTaskResult { candidates: alleles, evidence: ev }
}

/// The winning and near-winning complexes from the candidate complexes (ComplexCoverageCalculator and
/// ComplexCoverageRanking); `fas` the fragment alleles used, `per_gene` the product
#[allow(clippy::too_many_arguments)]
fn rank_product(config: &LilacConfig, rd: &RefData, product: &Product, fas: &[FragAlleles], recovered: &[AId], seqs: &[&SeqLoci]) -> Result<Vec<ComplexCoverage>> {
    let genes: Vec<HlaGene> = rd.gene_cache.gene_names.iter().copied().filter(|g| !g.is_pseudo()).collect();
    let totals: Vec<(u32, Totals)> = complexes::calculate_product(product, fas, &rd.alleles, &genes, threads(config), config.threads, config.top_score_threshold)?;
    let mut ranker = Ranker::new(rd, config.top_score_threshold, rd.genes.gene_count(), config.params.freq_score_penalty);
    // the complexes' alleles sorted by name
    let sorted_alleles = |i: u32| -> Vec<AId> {
        let mut v = product.alleles(i as usize);
        v.sort_by(|a, b| hmf_common::java_text::compare(rd.alleles.name(*a), rd.alleles.name(*b)));
        v
    };
    let mut items: Vec<Scored> = Vec::with_capacity(totals.len());
    for (i, t) in &totals {
        let al = sorted_alleles(*i);
        let mut s = ranker.score(&al, t.total(), recovered, seqs);
        s.index = *i;
        s.totals = *t;
        items.push(s);
    }
    let order = ranker.rank(&items, &|s: &Scored| sorted_alleles(s.index));
    // the kept complexes with their allele coverage as Java computes it
    let mut all: Vec<AId> = Vec::new();
    for cs in &product.per_gene {
        for c in cs {
            for a in c {
                if !all.contains(a) {
                    all.push(*a);
                }
            }
        }
    }
    let matrix = complexes::Matrix::new(fas, &all);
    Ok(order
        .into_iter()
        .map(|k| {
            let s = &items[k];
            let list = matrix.create(&product.alleles(s.index as usize));
            let mut c = ComplexCoverage::create(&rd.gene_cache, &rd.alleles, list);
            debug_assert_eq!((c.unique, c.shared, c.wild), (s.totals.unique, s.totals.shared, s.totals.wild));
            apply_score(&mut c, s);
            c
        })
        .collect())
}

/// LilacApplication.runGeneGroup
pub fn run_gene_group(config: &LilacConfig, rd: &mut RefData, ref_bases: Option<&RefBases>, source: Fragments, mut writers: Option<&mut GroupOut>) -> Result<GroupResult> {
    let started = std::time::Instant::now();
    let params = &config.params;
    let genes = rd.genes;
    let mut all_valid = true;
    let read_ctx = ReadCtx { params, ref_bases };
    let ref_bam_path = config.allele_calling_bam().to_string();
    crate::info!("finding read support in {} bam {}", if config.tumor_only() { "tumor" } else { "reference" }, ref_bam_path);

    let rd_ro: &RefData = rd;
    let factory = Factory::new(&rd_ro.gene_cache, &rd_ro.aa_with_inserts, &rd_ro.aa_with_deletes, ReadCtx { params, ref_bases });
    let mut ref_reader = GeneReader::new(&rd_ro.gene_cache, &params.hla_chr);
    let mut ref_bam: Option<Bam> = None;
    let (mut ref_frags, given_tumor, given_stop_loss) = match source {
        Fragments::Bams => {
            let mut b = Bam::open(&ref_bam_path, &config.ref_genome)?;
            let f = ref_reader.find_gene_fragments(&mut b, &factory, &read_ctx, rd_ro.stop_loss_indel.as_ref(), &rd_ro.pon_indels)?;
            ref_bam = Some(b);
            (f, None, Vec::new())
        }
        Fragments::Given { reference, tumor, stop_loss } => (reference, Some(tumor), stop_loss),
    };
    if !given_stop_loss.is_empty() {
        ref_reader.has_stop_loss = true;
        ref_reader.stop_loss_frags = given_stop_loss;
    }
    if let Some(e) = &rd_ro.gene_enrichment {
        for f in ref_frags.iter_mut() {
            check_additional_genes(e, &rd_ro.gene_cache, f);
        }
    }
    let median_base_qual: Option<i32> = if params.is_illumina() { Some(fragment::percentile_base_quality(ref_frags.iter(), BASE_QUAL_PERCENTILE) as i32) } else { None };
    if let Some(m) = median_base_qual {
        if m < LOW_BASE_QUAL_THRESHOLD as i32 {
            crate::info!("the median base quality ({m}) is below the low qual threshold ({LOW_BASE_QUAL_THRESHOLD})");
        }
    }
    let depth = fragment::gene_coverage(&rd_ro.gene_cache, ref_frags.iter());
    if !has_sufficient_depth(config, &depth) {
        write_failed(config, rd_ro, &depth, writers);
        return Ok(GroupResult::default());
    }
    if config.run_validation {
        let invalid = ref_frags.iter().filter(|f| !f.validate()).count();
        if invalid > 0 {
            crate::warn!("has {invalid} invalid fragments");
            all_valid = false;
        }
    }
    let mut pipeline = AaPipeline::new(&mut ref_frags);
    let ref_aa: Vec<usize> = pipeline.high_qual_ref.clone();
    crate::info!("totalFrags({})", ref_aa.len());

    // the gene tasks: their fragments are made in order (and each enriches the shared originals, LI3), then
    // the candidates and phasing run in parallel
    let mut raw = RawCounts::default();
    let mut task_frags: Vec<Vec<Fragment>> = Vec::new();
    for ctx in &rd_ro.contexts {
        task_frags.push(pipeline.reference_phasing_fragments(ctx, &rd_ro.gene_cache, &mut raw, params, !crate::compat::java("LI3")));
    }
    let raw_ref = &raw;
    let results: Vec<GeneTaskResult> = if threads(config) > 1 && rd_ro.contexts.len() > 1 {
        std::thread::scope(|s| {
            let handles: Vec<_> = rd_ro
                .contexts
                .iter()
                .zip(task_frags.iter())
                .map(|(ctx, frags)| s.spawn(move || gene_task(config, rd_ro, ctx, frags, raw_ref.aa.get(&ctx.gene))))
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        })
    } else {
        rd_ro.contexts.iter().zip(task_frags.iter()).map(|(ctx, frags)| gene_task(config, rd_ro, ctx, frags, raw_ref.aa.get(&ctx.gene))).collect()
    };
    drop(task_frags);
    let mut candidates: Vec<AId> = Vec::new();
    let mut combined_evidence: Vec<PhasedEvidence> = Vec::new();
    for r in &results {
        if r.candidates.is_empty() {
            continue;
        }
        if config.max_elimination_candidates == 0 || r.candidates.len() <= config.max_elimination_candidates as usize {
            candidates.extend(r.candidates.iter().copied());
        } else {
            let gene = rd_ro.alleles.gene(r.candidates[0]);
            // the cohort frequency alleles of the gene (Java adds uncached alleles here)
            for (name, _) in rd_ro.freqs.java_order() {
                if let Some(a) = rd_ro.alleles.find_four_digit(name) {
                    if rd_ro.alleles.gene(a) == gene {
                        candidates.push(a);
                    }
                }
            }
        }
    }
    for r in results {
        combined_evidence.extend(r.evidence);
    }

    let mut recovered: Vec<AId> = Vec::new();
    let mut stop_loss: Vec<(AId, Vec<u64>)> = Vec::new();
    if ref_reader.has_stop_loss {
        if let Some((_, a)) = rd_ro.stop_loss_alleles.first() {
            crate::info!("recovering stop loss allele({}) with {} fragments", rd_ro.alleles.name(*a), ref_reader.stop_loss_frags.len());
            stop_loss.push((*a, ref_reader.stop_loss_frags.clone()));
            if !candidates.contains(a) {
                recovered.push(*a);
            }
        }
    }
    if config.restricted_alleles.is_empty() {
        let mut missed: Vec<AId> = rd_ro.common_alleles.iter().copied().filter(|x| !candidates.contains(x) && !recovered.contains(x)).collect();
        if !missed.is_empty() {
            hmf_common::java_sort::java_sort(&mut missed, |a, b| rd_ro.alleles.compare(*a, *b));
            crate::info!("recovering {} common alleles: {}", missed.len(), rd_ro.alleles.join(&missed));
            recovered.extend(missed);
        }
    }
    candidates.extend(recovered.iter().copied());
    let missing_expected: Vec<AId> = rd_ro.actual_alleles.iter().copied().filter(|x| !candidates.contains(x)).collect();
    if !missing_expected.is_empty() {
        crate::info!("  including {} known actual alleles as candidates: {}", missing_expected.len(), rd_ro.alleles.join(&missing_expected));
        candidates.extend(missing_expected);
    }
    let cand_set: HashSet<AId> = candidates.iter().copied().collect();
    let mut candidate_seqs: Vec<&SeqLoci> = rd_ro.aa_seqs.iter().filter(|s| cand_set.contains(&s.allele)).collect();
    let ref_aa_counts = fragment::count_amino_acids(params, params.min_evidence_factor, ref_aa.iter().map(|i| &ref_frags[*i]));
    let ref_nuc_counts = fragment::count_nucleotides(params, params.min_evidence_factor, ref_aa.iter().map(|i| &ref_frags[*i]));
    let nuc_het = evidence::nucleotide_het_loci(&rd_ro.gene_cache, &ref_nuc_counts.heterozygous_loci());
    let mut candidate_nuc_seqs: Vec<&SeqLoci> = rd_ro.nuc_seqs.iter().filter(|s| cand_set.contains(&rd_ro.alleles.four_digit(s.allele))).collect();
    let recovered_set: HashSet<AId> = recovered.iter().copied().collect();
    let recovered_seqs: Vec<&SeqLoci> = rd_ro.aa_seqs.iter().filter(|s| recovered_set.contains(&s.allele)).collect();
    let aa_het = het_aa_loci(rd_ro, &pipeline.ref_aa_counts, &recovered_seqs);
    let mut mapper = Mapper { aa_het_loci: aa_het, nuc_het_loci: nuc_het, ref_nucleotides: pipeline.reference_nucleotides(), stop_loss, threads: threads(config) };
    let mut ref_fas = mapper.create_fragment_alleles(rd_ro, params, &mut ref_frags, &ref_aa, &candidate_seqs, &candidate_nuc_seqs);
    if ref_fas.is_empty() {
        crate::error!("failed to assign fragments to alleles");
        std::process::exit(1);
    }
    let mut hlay = rd_ro.hlay_nuc_seqs.as_ref().map(|s| HlaY::new(s.clone(), &mapper.aa_het_loci));
    if let Some(h) = hlay.as_mut() {
        h.check_threshold(config, &mut ref_fas, &mut ref_frags, &ref_aa);
    }
    let mut builder = ComplexBuilder::new(genes.gene_count());
    builder.filter_candidates(rd_ro, &ref_fas, &candidates, &recovered);
    if config.run_validation {
        let mut seen = HashSet::new();
        let dups = builder.unique_protein_alleles.iter().filter(|a| !seen.insert(rd_ro.alleles.name(**a))).count();
        if dups > 0 {
            crate::warn!("has {dups} duplicate alleles from complex building");
            all_valid = false;
        }
    }
    let confirmed_recovered = builder.confirmed_recovered_alleles.clone();
    let recovered_seqs: Vec<&SeqLoci> = recovered_seqs.into_iter().filter(|s| confirmed_recovered.contains(&s.allele)).collect();
    mapper.aa_het_loci = het_aa_loci(rd_ro, &pipeline.ref_aa_counts, &recovered_seqs);
    if let Some(h) = hlay.as_mut() {
        h.update_het_loci(&mapper.aa_het_loci);
    }
    let confirmed: HashSet<AId> = builder.unique_protein_alleles.iter().copied().collect();
    candidate_seqs.retain(|s| confirmed.contains(&s.allele));
    candidate_nuc_seqs.retain(|s| confirmed.contains(&rd_ro.alleles.four_digit(s.allele)));
    crate::debug!(
        "creating fragment alleles from amnioAcidfrags({}) candidateSequences({}) nucleotideSequences({})",
        ref_aa.len(),
        candidate_seqs.len(),
        candidate_nuc_seqs.len()
    );
    ref_fas = mapper.create_fragment_alleles(rd_ro, params, &mut ref_frags, &ref_aa, &candidate_seqs, &candidate_nuc_seqs);

    // down-sampling, then the wildcard filter of buildComplexes on the fragment alleles used (Java changes
    // the shared FragmentAlleles objects, and the reference list itself when it was not down-sampled)
    let calc_idx = downsample(config, &ref_fas, &ref_frags);
    let capped = calc_idx.len() < ref_fas.len();
    let calc_fas_before: Vec<FragAlleles> = calc_idx.iter().map(|i| ref_fas[*i].clone()).collect();
    let wildcards = coverage::find_wildcard_alleles(&rd_ro.alleles, &calc_fas_before);
    let unsupported = coverage::find_unsupported_wildcards(&rd_ro.alleles, &calc_fas_before, &wildcards);
    let mut calc_fas: Vec<FragAlleles> = Vec::with_capacity(calc_idx.len());
    let mut removed_idx: Vec<usize> = Vec::new();
    for &i in &calc_idx {
        ref_fas[i].remove_all(&unsupported);
        if ref_fas[i].full.is_empty() && ref_fas[i].wild.is_empty() {
            ref_frags[ref_fas[i].fi].set_scope(Scope::WildOnly, false);
            removed_idx.push(i);
        } else {
            calc_fas.push(ref_fas[i].clone());
        }
    }
    if !removed_idx.is_empty() {
        crate::debug!("removed frags({}) with wildcard alleles removed", removed_idx.len());
    }
    if !capped {
        let mut k = 0;
        ref_fas.retain(|_| {
            let keep = !removed_idx.contains(&k);
            k += 1;
            keep
        });
    }
    builder.unique_protein_alleles.retain(|a| !unsupported.contains(a));
    let per_gene: Vec<Vec<coverage::Complex>> = rd_ro
        .gene_cache
        .gene_names
        .iter()
        .filter(|g| !g.is_pseudo())
        .map(|g| ComplexBuilder::build_complexes_by_gene(&rd_ro.alleles, *g, &builder.unique_group_alleles, &builder.unique_protein_alleles))
        .collect();
    let simple: i64 = per_gene.iter().fold(1i64, |acc, x| acc.wrapping_mul(x.len() as i64));
    let product = if simple > COMPLEX_PERMS_THRESHOLD || simple < 0 {
        let common: Vec<AId> = rd_ro.common_alleles.iter().copied().filter(|a| builder.unique_protein_alleles.contains(a)).collect();
        let sizes: Vec<String> = rd_ro.gene_cache.gene_names.iter().filter(|g| !g.is_pseudo()).zip(per_gene.iter()).map(|(g, c)| format!("{}={}", g.long_name(), c.len())).collect();
        crate::info!("candidate permutations exceeds threshold, candidates({{{}}}) common({})", sizes.join(", "), common.len());
        let recovered_set: Vec<AId> = confirmed_recovered.clone();
        let mut top: Vec<AId> = Vec::new();
        for cs in &per_gene {
            top.extend(ranked_group_coverage(config, rd_ro, 10, &calc_fas, cs, &recovered_set));
        }
        for a in common {
            if !top.contains(&a) {
                top.push(a);
            }
        }
        let rejected: Vec<AId> = builder.unique_protein_alleles.iter().copied().filter(|a| !top.contains(a)).collect();
        crate::debug!("  discarding {} unlikely candidates: {}", rejected.len(), rd_ro.alleles.join(&rejected));
        Product { per_gene: ComplexBuilder::build_allele_complexes(rd_ro, &builder.unique_group_alleles, &top) }
    } else {
        Product { per_gene }
    };
    let ref_count_str = if ref_fas.len() > calc_fas.len() { format!("{} capped={}", ref_fas.len(), calc_fas.len()) } else { ref_fas.len().to_string() };
    crate::info!("calculating coverage for complexes({}) and ref alleles({ref_count_str})", product.len());
    let phase = std::time::Instant::now();
    let mut ranked = rank_product(config, rd_ro, &product, &calc_fas, &recovered, &candidate_seqs)?;
    crate::debug!("PerfStats: complex coverage and ranking {:.2}s", phase.elapsed().as_secs_f64());
    if ranked.is_empty() {
        crate::error!("failed to calculate complex coverage");
        std::process::exit(1);
    }
    if calc_fas.len() < ref_fas.len() {
        let filtered: Vec<coverage::Complex> = ranked.iter().map(|c| c.allele_ids()).collect();
        crate::debug!("recalculating coverage for complexes({}) and ref alleles({})", filtered.len(), ref_fas.len());
        let calculated = complexes::calculate_list(rd_ro, &filtered, &ref_fas, config.threads, config.top_score_threshold)?;
        let mut ranker = Ranker::new(rd_ro, 0.0, genes.gene_count(), params.freq_score_penalty);
        ranked = complexes::rank_coverages(&mut ranker, calculated, &recovered, &candidate_seqs);
    }
    // the winner is one of the ranked complexes (or the actual alleles' coverage) and is expanded in place, so
    // the candidates file shows it expanded, as in Java
    let mut win_pos: Option<usize> = Some(0);
    let mut forced: Option<ComplexCoverage> = None;
    if !rd_ro.actual_alleles.is_empty() {
        let actual = &rd_ro.actual_alleles;
        win_pos = ranked.iter().position(|c| c.alleles.len() == actual.len() && c.alleles.iter().all(|x| actual.contains(&x.allele)));
        if win_pos.is_none() {
            forced = Some(calc_protein_coverage(&rd_ro.gene_cache, &rd_ro.alleles, &ref_fas, actual));
        }
        crate::info!("forcing winning coverage to actual alleles");
    }
    let winning = match win_pos {
        Some(p) => {
            ranked[p].expand_to_six(&rd_ro.gene_cache, &rd_ro.alleles);
            ranked[p].clone()
        }
        None => {
            let mut w = forced.unwrap();
            w.expand_to_six(&rd_ro.gene_cache, &rd_ro.alleles);
            w
        }
    };
    let winning_alleles: Vec<AId> = winning.allele_ids();
    if let Some(h) = hlay.as_mut() {
        h.assign_reference(rd_ro, crate::log::debug_enabled(), config, &winning_alleles, &ref_fas, &ref_frags, &ref_aa);
    }
    let winning_seqs: Vec<&SeqLoci> = candidate_seqs.iter().copied().filter(|s| winning_alleles.contains(&s.allele)).collect();
    let winning_nuc_seqs: Vec<&SeqLoci> = candidate_nuc_seqs.iter().copied().filter(|s| winning_alleles.contains(&rd_ro.alleles.four_digit(s.allele))).collect();
    crate::info!("{}", info_header(genes));
    for c in &ranked {
        crate::info!("{}", c.as_line(&rd_ro.alleles, None));
    }
    let mut score_margin = 0.0;
    let mut next_info: Vec<String> = Vec::new();
    if ranked.len() > 1 {
        score_margin = ranked[0].score - ranked[1].score;
        for a in ranked[1].allele_ids() {
            if !winning_alleles.contains(&a) {
                next_info.push(rd_ro.alleles.name(a).to_string());
            }
        }
    }
    crate::info!("winning solution: {}", rd_ro.alleles.join(&winning_alleles));
    for fa in &ref_fas {
        let f = &mut ref_frags[fa.fi];
        if f.is_scope_set() {
            continue;
        }
        if winning_alleles.iter().any(|a| fa.contains(*a)) {
            f.set_scope(Scope::Solution, false);
        } else {
            f.set_scope(Scope::Candidate, false);
        }
    }
    let mut coding: Vec<SomaticCodingCount> = winning_alleles.iter().map(|a| SomaticCodingCount::new(*a)).collect();

    // extractTumorResults
    let mut tumor_cn: Vec<f64> = vec![0.0; 6];
    let mut somatic_count = 0;
    let tumor_cov: ComplexCoverage = if config.tumor_bam.is_empty() {
        ComplexCoverage::empty(&rd_ro.gene_cache)
    } else {
        let mut tumor_reader_state: Option<(GeneReader, Bam)> = None;
        let tc = if config.tumor_only() {
            winning.clone()
        } else {
            let mut tumor_frags = match &given_tumor {
                Some(t) => t.clone(),
                None => {
                    let mut tr = GeneReader::new(&rd_ro.gene_cache, &params.hla_chr);
                    let mut tb = Bam::open(&config.tumor_bam, &config.ref_genome)?;
                    let f = tr.find_gene_fragments(&mut tb, &factory, &read_ctx, rd_ro.stop_loss_indel.as_ref(), &rd_ro.pon_indels)?;
                    tumor_reader_state = Some((tr, tb));
                    f
                }
            };
            if let Some(e) = &rd_ro.gene_enrichment {
                for f in tumor_frags.iter_mut() {
                    check_additional_genes(e, &rd_ro.gene_cache, f);
                }
            }
            let n_tumor = tumor_frags.len();
            let tumor_idx = pipeline.comparison_coverage_fragments(&ref_frags, &mut tumor_frags, params);
            if let Some(w) = writers.as_deref_mut() {
                if config.writes(WriteType::Fragments) {
                    for i in &tumor_idx {
                        w.fragments.push(fragment_line("TUMOR", &tumor_frags[*i]));
                    }
                }
            }
            crate::info!("calculating tumor coverage from frags({} highQual={})", n_tumor, tumor_idx.len());
            let mut tumor_fas = mapper.create_fragment_alleles(rd_ro, params, &mut tumor_frags, &tumor_idx, &winning_seqs, &winning_nuc_seqs);
            if let Some(h) = hlay.as_mut() {
                if h.exceeds {
                    h.assign(rd_ro, crate::log::debug_enabled(), config, &winning_alleles, &mut tumor_fas, &tumor_frags, &tumor_idx, Source::Tumor);
                }
            }
            let mut tc = calc_protein_coverage(&rd_ro.gene_cache, &rd_ro.alleles, &tumor_fas, &winning_alleles);
            tc.populate_missing(&rd_ro.gene_cache, &winning_alleles);
            tc.expand_to_six(&rd_ro.gene_cache, &rd_ro.alleles);
            tc
        };
        if !config.copy_number_file.is_empty() {
            let cn = variant::load_copy_numbers(&config.copy_number_file, &rd_ro.gene_cache.gene_names);
            variant::assign_copy_numbers(rd_ro, &cn, &winning.alleles, &tc.alleles, &mut tumor_cn);
        }
        let annotation = VariantAnnotation::new(config, rd_ro)?;
        if !annotation.variants.is_empty() {
            somatic_count = annotation.variants.len() as i32;
            crate::info!("calculating somatic variant allele coverage");
            // the tumor reader: the reference reader itself in tumor-only mode (its filtered-read count goes
            // into the QC)
            let mut own: Option<(GeneReader, Bam)> = None;
            if !config.tumor_only() {
                own = Some(match tumor_reader_state.take() {
                    Some(x) => x,
                    None => (GeneReader::new(&rd_ro.gene_cache, &params.hla_chr), Bam::open(&config.tumor_bam, &config.ref_genome)?),
                });
            } else if ref_bam.is_none() {
                ref_bam = Some(Bam::open(&config.tumor_bam, &config.ref_genome)?);
            }
            for v in &annotation.variants {
                let mut frags = match own.as_mut() {
                    Some((r, b)) => variant::variant_fragments(r, b, &factory, &read_ctx, v)?,
                    None => variant::variant_fragments(&mut ref_reader, ref_bam.as_mut().unwrap(), &factory, &read_ctx, v)?,
                };
                let cov = annotation.assign_allele_coverage(rd_ro, config, v, &mut frags, &winning_seqs);
                let va: Vec<AId> = cov.iter().map(|x| x.allele).collect();
                crate::info!("  {} -> {}}}", v.display(), cov.iter().map(|x| x.display(&rd_ro.alleles)).collect::<Vec<_>>().join(", "));
                let indel = v.alt.len() != v.ref_.len();
                SomaticCodingCount::add_variant(&mut coding, &rd_ro.alleles, indel, v.effect, &va);
            }
        }
        tc
    };

    // extractRnaCoverage
    let rna_cov = if config.rna_bam.is_empty() {
        ComplexCoverage::empty(&rd_ro.gene_cache)
    } else {
        let mut rr = GeneReader::new(&rd_ro.gene_cache, &params.hla_chr);
        let mut rb = Bam::open(&config.rna_bam, &config.ref_genome)?;
        let mut rna_frags = rr.find_gene_fragments(&mut rb, &factory, &read_ctx, rd_ro.stop_loss_indel.as_ref(), &rd_ro.pon_indels)?;
        if let Some(e) = &rd_ro.gene_enrichment {
            for f in rna_frags.iter_mut() {
                check_additional_genes(e, &rd_ro.gene_cache, f);
            }
        }
        let n = rna_frags.len();
        let idx = pipeline.comparison_coverage_fragments(&ref_frags, &mut rna_frags, params);
        crate::info!("calculating RNA coverage from frags({} highQual={})", n, idx.len());
        let rna_fas = mapper.create_fragment_alleles(rd_ro, params, &mut rna_frags, &idx, &winning_seqs, &winning_nuc_seqs);
        let mut c = calc_protein_coverage(&rd_ro.gene_cache, &rd_ro.alleles, &rna_fas, &winning_alleles);
        c.populate_missing(&rd_ro.gene_cache, &winning_alleles);
        c.expand_to_six(&rd_ro.gene_cache, &rd_ro.alleles);
        c
    };

    if !all_valid {
        crate::error!("failed validation");
        write_failed(config, rd_ro, &depth, writers);
        return Ok(GroupResult::default());
    }
    let somatic_qc = SomaticVariantQC::create(somatic_count, &coding);
    let unmatched: Vec<&Fragment> = ref_aa.iter().map(|i| &ref_frags[*i]).filter(|f| f.scope.is_unmatched()).collect();
    let hlay_aa = rd_ro.hlay_aa_seqs.as_deref();
    let haplotype_qc = HaplotypeQC::create(params, &winning_seqs, hlay_aa, &combined_evidence, &ref_aa_counts, &unmatched);
    let aa_qc = AminoAcidQC::create(&winning_seqs, hlay_aa, &ref_aa_counts, &haplotype_qc.unmatched);
    let pon = ref_reader.unmatched_pon_indels(3);
    let unmatched_indels = ref_reader.unmatched_indels(3, rd_ro.stop_loss_indel.as_ref());
    for (k, v) in &unmatched_indels {
        crate::warn!("  UNMATCHED_INDEL - {v} fragments excluded with unmatched indel {k}");
    }
    for (k, v) in &pon {
        crate::debug!("  UNMATCHED_PON_INDEL - {v} fragments excluded with unmatched PON indel {k}");
    }
    let bam_qc = BamQC::new(&rd_ro.gene_cache, ref_reader.filtered_count, unmatched_indels.len() as i32, unmatched_indels.iter().map(|x| x.1).max().unwrap_or(0), &depth);
    let coverage_qc = CoverageQC::create(&rd_ro.gene_cache, &rd_ro.alleles, &ref_frags, &ref_aa, &winning);
    let hla_y_allele = hlay.as_ref().and_then(|h| h.selected_allele(&rd_ro.alleles)).map(|a| rd_ro.alleles.name(a).to_string());
    let summary = LilacQC::new(params, genes, score_margin, next_info.join(";"), median_base_qual, hla_y_allele, aa_qc, bam_qc, coverage_qc, haplotype_qc, somatic_qc);
    let ref_cov = if !config.tumor_only() { winning.clone() } else { ComplexCoverage::empty(&rd_ro.gene_cache) };
    let solution = SolutionSummary::create(&rd_ro.alleles, genes, ref_cov, tumor_cov, tumor_cn, &coding, rna_cov);

    // writeFileOutputs
    summary.log(&rd_ro.gene_cache, &config.sample);
    let summary_lines = solution.lines(&rd_ro.alleles);
    let qc_line = summary.body(&rd_ro.gene_cache).join("\t");
    if let Some(w) = writers.as_deref_mut() {
        if !config.output_dir.is_empty() {
            w.summary.extend(summary_lines.iter().cloned());
            w.qc.push(qc_line.clone());
            for c in &ranked {
                w.complexes.push(c.as_line(&rd_ro.alleles, Some(genes.name())));
            }
            if config.writes(WriteType::RefCounts) {
                w.files.push((config.form_file_id(LILAC_FILE_CANDIDATE_AA), ref_aa_counts.vertical_text(None)));
                w.files.push((config.form_file_id(LILAC_FILE_CANDIDATE_NUC), ref_nuc_counts.vertical_text(None)));
                for (g, c) in &pipeline.ref_aa_counts {
                    w.files.push((config.form_file_id(&format!("{}.aminoacids.txt", g.long_name())), c.vertical_text(raw.aa.get(g))));
                }
                for (g, c) in &pipeline.ref_nuc_counts {
                    w.files.push((config.form_file_id(&format!("{}.nucleotides.txt", g.long_name())), c.vertical_text(raw.nuc.get(g))));
                }
                if let Some(f) = hlay.as_ref().and_then(|h| h.allele_counts_file(config, &rd_ro.alleles)) {
                    w.files.push(f);
                }
            }
        }
        if config.writes(WriteType::Fragments) {
            for c in &ranked {
                let ids = c.allele_ids();
                let complex_str: Vec<&str> = ids.iter().map(|a| rd_ro.alleles.name(*a)).collect();
                for fa in coverage::filter_frag_alleles(&rd_ro.alleles, &ref_fas, &ids) {
                    let full: Vec<&str> = fa.full.iter().map(|a| rd_ro.alleles.name(*a)).collect();
                    let wild: Vec<&str> = fa.wild.iter().map(|a| rd_ro.alleles.name(*a)).collect();
                    let t = if fa.full.len() == 1 && fa.wild.is_empty() {
                        "UNIQUE_FULL"
                    } else if fa.full.is_empty() && fa.wild.len() == 1 {
                        "UNIQUE_WILD"
                    } else {
                        "SHARED"
                    };
                    let f = &ref_frags[fa.fi];
                    w.ref_fragments.push(format!("{}\t{}\t{}\t{}\t{}\t{}", complex_str.join(";"), f.id(), f.read_info(), t, full.join(";"), wild.join(";")));
                }
            }
            let src = if config.tumor_only() { "TUMOR" } else { "REFERENCE" };
            for f in &ref_frags {
                w.fragments.push(fragment_line(src, f));
            }
        }
        if let Some(h) = &hlay {
            if !h.fragment_lines.is_empty() && !config.output_dir.is_empty() {
                w.files.push((config.form_file_id(LILAC_FILE_HLA_Y_FRAGMENTS), crate::hlay::fragment_text(&h.fragment_lines)));
            }
        }
    }
    crate::info!("Lilac completed genes({}), mins({})", genes.name(), hmf_common::fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    let reference_alleles = solution.reference.as_ref().map(|c| c.allele_ids()).unwrap_or_default();
    Ok(GroupResult { ranked, summary_lines, qc_line: Some(qc_line), reference_alleles, qc: Some(summary) })
}

/// ComplexBuilder.rankedGroupCoverage: the top alleles of one gene's complexes
fn ranked_group_coverage(config: &LilacConfig, rd: &RefData, take: usize, fas: &[FragAlleles], complexes_g: &[coverage::Complex], recovered: &[AId]) -> Vec<AId> {
    let covs: Vec<ComplexCoverage> = complexes_g.iter().map(|c| calc_protein_coverage(&rd.gene_cache, &rd.alleles, fas, c)).collect();
    let mut ranker = Ranker::new(rd, 0.0, rd.genes.gene_count(), config.params.freq_score_penalty);
    let ranked = complexes::rank_coverages(&mut ranker, covs, recovered, &[]);
    let mut pairing: Vec<(AId, i32)> = Vec::new();
    let mut top: Vec<AId> = Vec::new();
    let count = |p: &Vec<(AId, i32)>, a: AId| p.iter().find(|x| x.0 == a).map(|x| x.1);
    for c in &ranked {
        let ids = c.allele_ids();
        let a1 = ids[0];
        if ids.len() == 1 {
            if !top.contains(&a1) {
                top.push(a1);
            }
        } else {
            let a2 = ids[1];
            let c1 = count(&pairing, a1);
            let c2 = count(&pairing, a2);
            if c1.is_some_and(|x| x >= 5) && c2.is_none() {
                continue;
            }
            if c2.is_some_and(|x| x >= 5) && c1.is_none() {
                continue;
            }
            match c1 {
                None => {
                    pairing.push((a1, 1));
                    if !top.contains(&a1) {
                        top.push(a1);
                    }
                }
                Some(x) => pairing.iter_mut().find(|p| p.0 == a1).unwrap().1 = x + 1,
            }
            match c2 {
                None => {
                    pairing.push((a2, 1));
                    if !top.contains(&a2) {
                        top.push(a2);
                    }
                }
                Some(x) => pairing.iter_mut().find(|p| p.0 == a2).unwrap().1 = x + 1,
            }
        }
        if top.len() >= take {
            break;
        }
    }
    top
}

/// HlaComplexFile.infoHeader
fn info_header(genes: GeneSelector) -> String {
    let mut v: Vec<String> = [
        "Score",
        "ComplexityPenalty",
        "Complexity",
        "HomozygousCount",
        "CohortFrequencyPenalty",
        "CohortFrequency",
        "RecoveryCount",
        "WildcardCount",
        "TotalCoverage",
        "UniqueCoverage",
        "SharedCoverage",
        "WildCoverage",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for g in genes.genes() {
        v.push(format!("{}_1", g.short_name()));
        v.push(format!("{}_2", g.short_name()));
    }
    v.join("\t")
}
