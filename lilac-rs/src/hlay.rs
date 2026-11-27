//! HlaYCoverage: fragments that match the HLA-Y pseudogene.

use hmf_common::java_hash::{bucket_index, grown_capacity};

use crate::config::LilacConfig;
use crate::coverage::FragAlleles;
use crate::fragment::{Fragment, Scope};
use crate::hla::{AId, Alleles, HlaGene};
use crate::refdata::RefData;
use crate::seq::{MatchType, SeqLoci};

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Source {
    Reference,
    Tumor,
    Rna,
}

impl Source {
    pub fn name(self) -> &'static str {
        match self {
            Source::Reference => "REFERENCE",
            Source::Tumor => "TUMOR",
            Source::Rna => "RNA",
        }
    }
}

const UNIQUE: usize = 0;
const SHARED_HLAY: usize = 1;
const SHARED: usize = 2;

pub struct HlaY {
    seqs: Vec<SeqLoci>,
    het_loci: Vec<i32>,
    removed_ref: Vec<FragAlleles>,
    pub exceeds: bool,
    /// per source (insertion order: identity-hashed keys under -XX:hashCode=2), the HLA-Y alleles in
    /// first-insertion order with their counts
    counts: Vec<(Source, Vec<(AId, [i32; 3])>)>,
    misc: Vec<(Source, [i32; 2])>,
    /// HLA_Y.fragments.tsv lines (written in debug mode only)
    pub fragment_lines: Vec<String>,
}

fn matches_full(s: &SeqLoci, f: &Fragment) -> bool {
    let loci: Vec<i32> = f.nucs.keys().collect();
    let nucs = f.nucleotides_str(loci.iter().copied());
    s.match_type_str(&nucs, &loci) == MatchType::Full
}

impl HlaY {
    pub fn new(seqs: Vec<SeqLoci>, aa_het: &[(HlaGene, std::collections::BTreeMap<i32, Vec<crate::tok::Tok>>)]) -> HlaY {
        let het_loci = aa_het.iter().find(|x| x.0 == HlaGene::A).map(|x| x.1.keys().copied().collect()).unwrap_or_default();
        HlaY { seqs, het_loci, removed_ref: Vec::new(), exceeds: false, counts: Vec::new(), misc: Vec::new(), fragment_lines: Vec::new() }
    }

    pub fn update_het_loci(&mut self, aa_het: &[(HlaGene, std::collections::BTreeMap<i32, Vec<crate::tok::Tok>>)]) {
        self.het_loci = aa_het.iter().find(|x| x.0 == HlaGene::A).map(|x| x.1.keys().copied().collect()).unwrap_or_default();
    }

    fn has_het_locus(&self, f: &Fragment) -> bool {
        f.aas.keys().any(|l| self.het_loci.binary_search(&l).is_ok())
    }

    /// checkThreshold: takes the HLA-Y fragments out of the reference fragment alleles when there are enough
    pub fn check_threshold(&mut self, config: &LilacConfig, fas: &mut Vec<FragAlleles>, frags: &mut [Fragment], subset: &[usize]) {
        let mut unique = 0;
        let mut matched: Vec<usize> = Vec::new(); // positions in fas
        for &fi in subset {
            if !self.has_het_locus(&frags[fi]) {
                continue;
            }
            let mut found = false;
            let mut mf: Option<usize> = None;
            for s in &self.seqs {
                if matches_full(s, &frags[fi]) {
                    found = true;
                    let id = frags[fi].id();
                    mf = fas.iter().position(|x| frags[x.fi].id() == id);
                    break;
                }
            }
            if !found {
                continue;
            }
            match mf {
                None => {
                    unique += 1;
                    frags[fi].set_scope(Scope::HlaY, true);
                }
                Some(p) => matched.push(p),
            }
        }
        let total = unique + matched.len();
        let threshold = subset.len() as f64 * config.hla_y_percent_threshold;
        self.exceeds = unique as f64 >= threshold;
        if total > 0 {
            crate::info!(
                "HLA-Y fragments({} unique={}) shared={}) {}",
                total,
                unique,
                matched.len(),
                if self.exceeds { "above threshold" } else { "below threshold" }
            );
            if self.exceeds {
                // matchedFragmentAlleles.forEach(fragAlleles::remove): each object once, first occurrence
                let objs: Vec<FragAlleles> = matched.iter().map(|p| fas[*p].clone()).collect();
                let mut remove: Vec<usize> = matched.clone();
                remove.sort();
                remove.dedup();
                for p in remove.into_iter().rev() {
                    fas.remove(p);
                }
                self.removed_ref.extend(objs);
                for fa in &self.removed_ref {
                    frags[fa.fi].set_scope(Scope::HlaY, true);
                }
            }
        }
    }

    pub fn assign_reference(&mut self, rd: &RefData, debug: bool, config: &LilacConfig, winners: &[AId], fas: &[FragAlleles], frags: &[Fragment], subset: &[usize]) {
        let mut all: Vec<FragAlleles> = fas.to_vec();
        all.extend(self.removed_ref.iter().cloned());
        self.assign(rd, debug, config, winners, &mut all, frags, subset, Source::Reference);
    }

    /// assignFragments
    #[allow(clippy::too_many_arguments)]
    pub fn assign(&mut self, rd: &RefData, debug: bool, config: &LilacConfig, winners: &[AId], fas: &mut Vec<FragAlleles>, frags: &[Fragment], subset: &[usize], source: Source) {
        let mut counts: Vec<(AId, [i32; 3])> = Vec::new();
        let mut misc = [0i32; 2];
        let mut matched: Vec<usize> = Vec::new();
        let mut unique = 0;
        let exon = rd.gene_cache.aa_boundaries(HlaGene::A);
        for &fi in subset {
            let f = &frags[fi];
            let aa_loci: Vec<i32> = f.aas.keys().filter(|l| self.het_loci.binary_search(l).is_ok()).collect();
            if aa_loci.is_empty() {
                continue;
            }
            let alleles: Vec<AId> = self.seqs.iter().filter(|s| matches_full(s, f)).map(|s| s.allele).collect();
            if alleles.is_empty() {
                continue;
            }
            // updateMiscCounts
            if f.aas.contains(227) && f.aa(227) == b'X' as crate::tok::Tok {
                misc[0] += 1;
            }
            let (s3, e3) = (exon[1] + 1, exon[2]);
            if f.min_aa_locus() >= s3 && f.max_aa_locus() <= e3 {
                misc[1] += 1;
            }
            let id = f.id();
            let m = fas.iter().position(|x| frags[x.fi].id() == id && winners.iter().any(|a| x.contains(*a)));
            match m {
                Some(p) => matched.push(p),
                None => unique += 1,
            }
            if debug && self.exceeds && !config.output_dir.is_empty() {
                let names: Vec<&str> = alleles.iter().map(|a| rd.alleles.name(*a)).collect();
                let loci: Vec<String> = aa_loci.iter().map(|x| x.to_string()).collect();
                self.fragment_lines.push(format!("{}\t{}\t{}\t{}\t{}", id, f.read_info(), names.join(";"), m.is_some(), loci.join(";")));
            }
            let multi = alleles.len() > 1;
            for a in alleles {
                let e = match counts.iter_mut().find(|x| x.0 == a) {
                    Some(e) => e,
                    None => {
                        counts.push((a, [0; 3]));
                        counts.last_mut().unwrap()
                    }
                };
                if m.is_some() {
                    e.1[SHARED] += 1;
                } else if multi {
                    e.1[SHARED_HLAY] += 1;
                } else {
                    e.1[UNIQUE] += 1;
                }
            }
        }
        match self.counts.iter_mut().find(|x| x.0 == source) {
            Some(x) => x.1 = counts,
            None => self.counts.push((source, counts)),
        }
        match self.misc.iter_mut().find(|x| x.0 == source) {
            Some(x) => x.1 = misc,
            None => self.misc.push((source, misc)),
        }
        if source != Source::Reference {
            crate::info!("HLA-Y src({}) fragments({} unique={}) shared={})", source.name(), unique + matched.len(), unique, matched.len());
            if self.exceeds {
                let mut remove = matched.clone();
                remove.sort();
                remove.dedup();
                for p in remove.into_iter().rev() {
                    fas.remove(p);
                }
            }
        }
    }

    fn java_order(alleles: &Alleles, v: &[(AId, [i32; 3])]) -> Vec<(AId, [i32; 3])> {
        let cap = grown_capacity(16, v.len());
        let mut k: Vec<(usize, (AId, [i32; 3]))> = v.iter().map(|x| (bucket_index(alleles.hash(x.0), cap), *x)).collect();
        k.sort_by_key(|x| x.0);
        k.into_iter().map(|x| x.1).collect()
    }

    /// getSelectedAllele: the HLA-Y allele with the most unique reference fragments
    pub fn selected_allele(&self, alleles: &Alleles) -> Option<AId> {
        if !self.exceeds {
            return None;
        }
        let r = self.counts.iter().find(|x| x.0 == Source::Reference)?;
        let mut best: Option<AId> = None;
        let mut max = 0;
        for (a, c) in Self::java_order(alleles, &r.1) {
            if c[UNIQUE] > 0 && c[UNIQUE] > max {
                max = c[UNIQUE];
                best = Some(a);
            }
        }
        best
    }

    /// writeAlleleCounts: HLA_Y.coverage.tsv (path and content), when above the threshold
    pub fn allele_counts_file(&self, config: &LilacConfig, alleles: &Alleles) -> Option<(String, String)> {
        if config.output_dir.is_empty() || !self.exceeds {
            return None;
        }
        for (s, m) in &self.misc {
            let outcome = if m[0] == 0 && m[1] > 0 { "NOVEL" } else { "REF" };
            crate::debug!("HLA-Y_MISC_COUNTS:{},{},{},{},{}", config.sample, s.name(), m[0], m[1], outcome);
        }
        let mut text = "Source\tAllele\tTotal\tShared\tHlaYShared\tUnique\n".to_string();
        for (s, v) in &self.counts {
            for (a, c) in Self::java_order(alleles, v) {
                text.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", s.name(), alleles.name(a), c[0] + c[1] + c[2], c[SHARED], c[SHARED_HLAY], c[UNIQUE]));
            }
        }
        Some((config.form_file_id(crate::constants::LILAC_FILE_HLA_Y_COVERAGE), text))
    }
}

/// The HLA-Y fragments file: Java never closes it, so only whole 8192-character buffers reach the disk (LI13)
pub fn fragment_text(lines: &[String]) -> String {
    let mut text = String::from("ReadId\tReadInfo\tAlleles\tShared\tLoci\n");
    for l in lines {
        text.push_str(l);
        text.push('\n');
    }
    if crate::compat::java("LI13") {
        let keep = text.len() / 8192 * 8192;
        text.truncate(keep);
    }
    text
}
