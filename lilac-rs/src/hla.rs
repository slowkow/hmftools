//! HlaGene, GeneSelector, HlaAllele and HlaAlleleCache.
//!
//! An allele is an index into the `Alleles` table of the gene group being solved. The Java cache hands out
//! one object per group and per four-digit (protein) allele, and the code compares these by identity in a
//! few places; here each cached allele has one index, so index equality is that identity. Alleles requested
//! at six or eight digits (the nucleotide sequences) get a new entry on every request, as Java creates a new
//! object each time; they are only used through their gene and their four-digit allele.

use std::cmp::Ordering;
use std::collections::HashMap;

use hmf_common::java_hash::java_string_hash;

#[derive(Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum HlaGene {
    A,
    B,
    C,
    Y,
    H,
    Dqb1,
    Dpa1,
    Dpb1,
    Dqa1,
    Drb1,
    Drb3,
    Drb4,
    Drb5,
    None,
}

impl HlaGene {
    pub const ALL: [HlaGene; 14] = [
        HlaGene::A,
        HlaGene::B,
        HlaGene::C,
        HlaGene::Y,
        HlaGene::H,
        HlaGene::Dqb1,
        HlaGene::Dpa1,
        HlaGene::Dpb1,
        HlaGene::Dqa1,
        HlaGene::Drb1,
        HlaGene::Drb3,
        HlaGene::Drb4,
        HlaGene::Drb5,
        HlaGene::None,
    ];

    pub fn ordinal(self) -> usize {
        self as usize
    }
    pub fn is_class1(self) -> bool {
        matches!(self, HlaGene::A | HlaGene::B | HlaGene::C | HlaGene::Y | HlaGene::H | HlaGene::None)
    }
    pub fn is_pseudo(self) -> bool {
        matches!(self, HlaGene::Y | HlaGene::H)
    }
    pub fn is_debug(self) -> bool {
        self == HlaGene::None
    }
    pub fn has_frequencies(self) -> bool {
        !matches!(self, HlaGene::Drb3 | HlaGene::Drb4 | HlaGene::Drb5)
    }
    /// the enum constant name, e.g. HLA_A
    pub fn enum_name(self) -> &'static str {
        match self {
            HlaGene::A => "HLA_A",
            HlaGene::B => "HLA_B",
            HlaGene::C => "HLA_C",
            HlaGene::Y => "HLA_Y",
            HlaGene::H => "HLA_H",
            HlaGene::Dqb1 => "HLA_DQB1",
            HlaGene::Dpa1 => "HLA_DPA1",
            HlaGene::Dpb1 => "HLA_DPB1",
            HlaGene::Dqa1 => "HLA_DQA1",
            HlaGene::Drb1 => "HLA_DRB1",
            HlaGene::Drb3 => "HLA_DRB3",
            HlaGene::Drb4 => "HLA_DRB4",
            HlaGene::Drb5 => "HLA_DRB5",
            HlaGene::None => "NONE",
        }
    }
    /// toString / longName: HLA-A, "" for NONE
    pub fn long_name(self) -> String {
        if self == HlaGene::None {
            return String::new();
        }
        self.enum_name().replace('_', "-")
    }
    /// shortName: A, DQB1; "" for NONE
    pub fn short_name(self) -> &'static str {
        if self == HlaGene::None {
            return "";
        }
        &self.enum_name()[4..]
    }
    /// HlaGene.fromString: "HLA-A" or "A"; None for anything else
    pub fn from_string(s: &str) -> Option<HlaGene> {
        let g = s.strip_prefix("HLA-").unwrap_or(s);
        let name = format!("HLA_{g}");
        HlaGene::ALL.iter().copied().find(|x| x.enum_name() == name)
    }
}

/// GeneSelector: the gene groups solved one after the other
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum GeneSelector {
    MhcClass1,
    HlaDqb1,
    HlaDpa1,
    HlaDpb1,
    HlaDqa1,
    HlaDrb1,
    HlaDrb,
}

impl GeneSelector {
    pub const ALL: [GeneSelector; 7] = [
        GeneSelector::MhcClass1,
        GeneSelector::HlaDqb1,
        GeneSelector::HlaDpa1,
        GeneSelector::HlaDpb1,
        GeneSelector::HlaDqa1,
        GeneSelector::HlaDrb1,
        GeneSelector::HlaDrb,
    ];

    pub fn name(self) -> &'static str {
        match self {
            GeneSelector::MhcClass1 => "MHC_CLASS_1",
            GeneSelector::HlaDqb1 => "HLA_DQB1",
            GeneSelector::HlaDpa1 => "HLA_DPA1",
            GeneSelector::HlaDpb1 => "HLA_DPB1",
            GeneSelector::HlaDqa1 => "HLA_DQA1",
            GeneSelector::HlaDrb1 => "HLA_DRB1",
            GeneSelector::HlaDrb => "HLA_DRB",
        }
    }
    pub fn included_in_all(self) -> bool {
        !matches!(self, GeneSelector::HlaDrb1 | GeneSelector::HlaDrb)
    }
    /// the genes solved
    pub fn genes(self) -> &'static [HlaGene] {
        match self {
            GeneSelector::MhcClass1 => &[HlaGene::A, HlaGene::B, HlaGene::C],
            GeneSelector::HlaDqb1 => &[HlaGene::Dqb1],
            GeneSelector::HlaDpa1 => &[HlaGene::Dpa1],
            GeneSelector::HlaDpb1 => &[HlaGene::Dpb1],
            GeneSelector::HlaDqa1 => &[HlaGene::Dqa1],
            GeneSelector::HlaDrb1 => &[HlaGene::Drb1],
            GeneSelector::HlaDrb => &[HlaGene::Drb1, HlaGene::Drb3, HlaGene::Drb4, HlaGene::Drb5],
        }
    }
    /// the genes solved and their pseudogenes
    pub fn all_genes(self) -> &'static [HlaGene] {
        match self {
            GeneSelector::MhcClass1 => &[HlaGene::A, HlaGene::B, HlaGene::C, HlaGene::Y, HlaGene::H],
            _ => self.genes(),
        }
    }
    pub fn gene_count(self) -> usize {
        self.genes().len()
    }
    pub fn contains(self, gene: HlaGene) -> bool {
        self.genes().contains(&gene)
    }
    pub fn covers_mhc_class1(self) -> bool {
        self.contains(HlaGene::A)
    }
    pub fn value_of(s: &str) -> Option<GeneSelector> {
        GeneSelector::ALL.iter().copied().find(|g| g.name() == s)
    }

    /// GeneSelector.parseArg. Java collects the groups in a HashSet of enum constants, so the order in which
    /// they are solved (and written) follows their identity hashes, which change from run to run (LI1); this
    /// gives the order of a JVM with constant identity hashes (-XX:hashCode=2): first-mention order.
    pub fn parse_arg(arg: &str) -> Result<Vec<GeneSelector>, String> {
        let mut genes: Vec<GeneSelector> = Vec::new();
        let add = |g: GeneSelector, genes: &mut Vec<GeneSelector>| {
            if !genes.contains(&g) {
                genes.push(g);
            }
        };
        for group in hmf_common::java_text::split(arg, ',') {
            if group.eq_ignore_ascii_case("ALL") {
                for g in GeneSelector::ALL {
                    if g.included_in_all() {
                        add(g, &mut genes);
                    }
                }
                continue;
            }
            match GeneSelector::value_of(group) {
                Some(g) => add(g, &mut genes),
                None => {
                    return Err(format!(
                        "java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.lilac.GeneSelector.{group}"
                    ));
                }
            }
        }
        Ok(genes)
    }
}

pub type AId = u32;

#[derive(Clone, Debug)]
pub struct AlleleInfo {
    pub gene: HlaGene,
    pub group: String,
    pub protein: String,
    pub synonymous: String,
    pub synonymous_non_coding: String,
    pub name: String,
    pub hash: i32,
    pub protein_number: i32,
    /// mFourDigit (None: this allele)
    pub four_digit: Option<AId>,
    /// mGroup (None: this allele)
    pub group_allele: Option<AId>,
    pub has_wildcards: bool,
}

/// The parts of an allele name, as HlaAllele.fromString splits it
pub struct AlleleParts {
    pub gene: HlaGene,
    pub group: String,
    pub protein: String,
    pub synonymous: String,
    pub synonymous_non_coding: String,
}

pub fn parse_allele(s: &str) -> Option<AlleleParts> {
    let star = s.find('*')?;
    let gene = HlaGene::from_string(&s[..star])?;
    let rest = &s[star + 1..];
    let split = hmf_common::java_text::split(rest, ':');
    let part = |i: usize| -> String { split.get(i).map(|x| x.to_string()).unwrap_or_default() };
    Some(AlleleParts { gene, group: part(0), protein: part(1), synonymous: part(2), synonymous_non_coding: part(3) })
}

pub fn allele_name(gene: HlaGene, group: &str, protein: &str, syn: &str, syn_nc: &str) -> String {
    let g = gene.short_name();
    if protein.is_empty() {
        format!("{g}*{group}")
    } else if syn.is_empty() {
        format!("{g}*{group}:{protein}")
    } else if syn_nc.is_empty() {
        format!("{g}*{group}:{protein}:{syn}")
    } else {
        format!("{g}*{group}:{protein}:{syn}:{syn_nc}")
    }
}

fn protein_number(protein: &str) -> i32 {
    if protein.is_empty() {
        return 0;
    }
    match hmf_common::java_text::parse_int(protein) {
        Some(v) => v,
        None => hmf_common::java_text::parse_int(&protein[..protein.len() - 1]).unwrap_or(0),
    }
}

/// HlaAlleleCache, plus the allele objects it hands out
#[derive(Clone, Default)]
pub struct Alleles {
    pub v: Vec<AlleleInfo>,
    groups: HashMap<String, AId>,
    four_digits: HashMap<String, AId>,
}

impl Alleles {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn get(&self, a: AId) -> &AlleleInfo {
        &self.v[a as usize]
    }
    #[inline]
    pub fn name(&self, a: AId) -> &str {
        &self.v[a as usize].name
    }
    #[inline]
    pub fn gene(&self, a: AId) -> HlaGene {
        self.v[a as usize].gene
    }
    #[inline]
    pub fn hash(&self, a: AId) -> i32 {
        self.v[a as usize].hash
    }
    #[inline]
    pub fn four_digit(&self, a: AId) -> AId {
        self.v[a as usize].four_digit.unwrap_or(a)
    }
    #[inline]
    pub fn group(&self, a: AId) -> AId {
        self.v[a as usize].group_allele.unwrap_or(a)
    }
    #[inline]
    pub fn has_wildcards(&self, a: AId) -> bool {
        self.v[a as usize].has_wildcards
    }
    pub fn set_has_wildcards(&mut self, a: AId) {
        self.v[a as usize].has_wildcards = true;
    }

    fn push(&mut self, gene: HlaGene, group: &str, protein: &str, syn: &str, syn_nc: &str, four: Option<AId>, grp: Option<AId>) -> AId {
        let name = allele_name(gene, group, protein, syn, syn_nc);
        let hash = java_string_hash(&name);
        self.v.push(AlleleInfo {
            gene,
            group: group.to_string(),
            protein: protein.to_string(),
            synonymous: syn.to_string(),
            synonymous_non_coding: syn_nc.to_string(),
            name,
            hash,
            protein_number: protein_number(protein),
            four_digit: four,
            group_allele: grp,
            has_wildcards: false,
        });
        (self.v.len() - 1) as AId
    }

    fn get_or_create_group(&mut self, gene: HlaGene, group: &str) -> AId {
        let name = allele_name(gene, group, "", "", "");
        if let Some(g) = self.groups.get(&name) {
            return *g;
        }
        let id = self.push(gene, group, "", "", "", None, None);
        self.groups.insert(name, id);
        id
    }

    fn get_or_create_four_digit(&mut self, gene: HlaGene, group: &str, protein: &str, group_id: AId) -> AId {
        let name = allele_name(gene, group, protein, "", "");
        if let Some(f) = self.four_digits.get(&name) {
            return *f;
        }
        let id = self.push(gene, group, protein, "", "", None, Some(group_id));
        self.four_digits.insert(name, id);
        id
    }

    /// HlaAlleleCache.request: a new allele object (with the cached four-digit and group alleles) every time
    pub fn request(&mut self, s: &str) -> Option<AId> {
        let p = parse_allele(s)?;
        let g = self.get_or_create_group(p.gene, &p.group);
        // tmp.asFourDigit(): the allele itself when it has no synonymous field
        let f = self.get_or_create_four_digit(p.gene, &p.group, &p.protein, g);
        Some(self.push(p.gene, &p.group, &p.protein, &p.synonymous, &p.synonymous_non_coding, Some(f), Some(g)))
    }

    /// HlaAlleleCache.requestFourDigit: the cached four-digit allele
    pub fn request_four_digit(&mut self, s: &str) -> Option<AId> {
        let p = parse_allele(s)?;
        let g = self.get_or_create_group(p.gene, &p.group);
        Some(self.get_or_create_four_digit(p.gene, &p.group, &p.protein, g))
    }

    /// findFourDigitAllele: an existing four-digit allele by name
    pub fn find_four_digit(&self, s: &str) -> Option<AId> {
        self.four_digits.get(s).copied()
    }

    /// HlaAllele.compareTo. Java returns "equal" for any two alleles with the same String hash code, and
    /// some allele names share one (LI14); the correction compares the names' fields then.
    pub fn compare(&self, a: AId, b: AId) -> Ordering {
        let x = self.get(a);
        let y = self.get(b);
        if x.hash == y.hash && (x.name == y.name || crate::compat::java("LI14")) {
            return Ordering::Equal;
        }
        let c = x.gene.ordinal().cmp(&y.gene.ordinal());
        if c != Ordering::Equal {
            return c;
        }
        let c = hmf_common::java_text::compare(&x.group, &y.group);
        if c != Ordering::Equal {
            return c;
        }
        if x.protein_number != y.protein_number {
            return if x.protein_number > y.protein_number { Ordering::Greater } else { Ordering::Less };
        }
        let c = hmf_common::java_text::compare(&x.synonymous, &y.synonymous);
        if c != Ordering::Equal {
            return c;
        }
        hmf_common::java_text::compare(&x.synonymous_non_coding, &y.synonymous_non_coding)
    }

    /// HlaAllele.toString(list): names joined with ", "
    pub fn join(&self, alleles: &[AId]) -> String {
        alleles.iter().map(|a| self.name(*a)).collect::<Vec<_>>().join(", ")
    }
    pub fn join_max(&self, alleles: &[AId], max: usize) -> String {
        let n = if max > 0 { alleles.len().min(max) } else { alleles.len() };
        self.join(&alleles[..n])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genes() {
        assert_eq!(HlaGene::from_string("HLA-A"), Some(HlaGene::A));
        assert_eq!(HlaGene::from_string("DQB1"), Some(HlaGene::Dqb1));
        assert_eq!(HlaGene::from_string("DQB1-AS1"), None);
        assert_eq!(HlaGene::from_string("NONE"), None);
        assert_eq!(HlaGene::Dpb1.long_name(), "HLA-DPB1");
        assert_eq!(HlaGene::Dpb1.short_name(), "DPB1");
        assert_eq!(GeneSelector::parse_arg("ALL").unwrap().len(), 5);
        assert_eq!(GeneSelector::parse_arg("HLA_DQB1,MHC_CLASS_1,HLA_DQB1").unwrap(), vec![GeneSelector::HlaDqb1, GeneSelector::MhcClass1]);
        assert!(GeneSelector::parse_arg("mhc_class_1").is_err());
    }

    #[test]
    fn alleles() {
        let mut a = Alleles::new();
        let x = a.request("A*01:01:01:02").unwrap();
        let y = a.request("A*01:01:02").unwrap();
        let f = a.request_four_digit("A*01:01").unwrap();
        assert_eq!(a.four_digit(x), f);
        assert_eq!(a.four_digit(y), f);
        assert_eq!(a.name(x), "A*01:01:01:02");
        assert_eq!(a.name(a.group(f)), "A*01");
        let b = a.request_four_digit("A*01:100").unwrap();
        let c = a.request_four_digit("A*01:20").unwrap();
        // protein numbers compare as numbers, names as strings
        assert_eq!(a.compare(c, b), Ordering::Less);
        // two names with one String hash code: equal in Java (LI14), ordered by protein number when corrected
        let p = a.request_four_digit("A*32:11Q").unwrap();
        let q = a.request_four_digit("A*32:122").unwrap();
        assert_eq!(a.hash(p), a.hash(q));
        assert_eq!(crate::compat::with_java(true, || a.compare(p, q)), Ordering::Equal);
        assert_eq!(crate::compat::with_java(false, || a.compare(p, q)), Ordering::Less);
        let n = a.request_four_digit("C*04:09N").unwrap();
        assert_eq!(a.get(n).protein_number, 9);
    }
}
