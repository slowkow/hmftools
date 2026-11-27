//! hmf-common SmallVariantFactory: the small variants of a PURPLE (or PAVE) VCF as SmallVariant objects,
//! with what VariantContextDecorator reads from each record (the PAVE impact, PURPLE's copy numbers,
//! hotspot, ClinVar, phase sets, reportable transcripts), VariantTranscriptImpact (PAVE_TI),
//! AltTranscriptReportableInfo and VariantEffect. ORANGE reads its somatic and germline variants this way.
//!
//! As Java:
//! - the factory's filters: PASS only for `pass_only` (htsjdk PassingVariantFilter: unfiltered or PASS),
//!   always HumanChromosomeFilter (1-22, X and Y, with or without "chr") and NTFilter (no NT attribute, or
//!   NT=ref);
//! - a record is kept when the tumor's genotype has an AD of two or more values and a total read count
//!   above zero (the larger of DP and the AD sum);
//! - the sample checks: the tumor, a reference (when given) and an RNA sample (when given) must be in
//!   the header, and the header must declare the AD format field;
//! - the INFO accessors are htsjdk's getAttributeAs* (`variant_context::Attributes`), which throw on a
//!   value that does not parse: here a panic, as for every port reading a VCF this way.
//!
//! Enum-valued fields are kept as their names.

use anyhow::{anyhow, bail, Result};

use crate::variant_context::{Genotype, GenotypeType, VcfHeader, VcfReader, VcfRecord, VariantType};

/// CodingEffect names, in the enum's order
pub const CODING_EFFECTS: [&str; 6] = ["NONSENSE_OR_FRAMESHIFT", "SPLICE", "MISSENSE", "SYNONYMOUS", "NONE", "UNDEFINED"];

/// VariantEffect: (enum name, effect string), in the enum's order
pub const VARIANT_EFFECTS: [(&str, &str); 20] = [
    ("STOP_GAINED", "stop_gained"),
    ("STOP_LOST", "stop_lost"),
    ("START_LOST", "start_lost"),
    ("FRAMESHIFT", "frameshift_variant"),
    ("SPLICE_ACCEPTOR", "splice_acceptor_variant"),
    ("SPLICE_DONOR", "splice_donor_variant"),
    ("INFRAME_INSERTION", "inframe_insertion"),
    ("INFRAME_DELETION", "inframe_deletion"),
    ("MISSENSE", "missense_variant"),
    ("PHASED_INFRAME_INSERTION", "phased_inframe_insertion"),
    ("PHASED_INFRAME_DELETION", "phased_inframe_deletion"),
    ("PHASED_MISSENSE", "phased_missense"),
    ("SYNONYMOUS", "synonymous_variant"),
    ("PHASED_SYNONYMOUS", "phased_synonymous"),
    ("INTRONIC", "intron_variant"),
    ("FIVE_PRIME_UTR", "5_prime_UTR_variant"),
    ("THREE_PRIME_UTR", "3_prime_UTR_variant"),
    ("UPSTREAM_GENE", "upstream_gene_variant"),
    ("NON_CODING_TRANSCRIPT", "non_coding_transcript_exon_variant"),
    ("OTHER", "other"),
];

/// VariantEffect.effectsToList: the enum names of the '&'-separated effect strings (an unknown one is OTHER)
pub fn effects_to_list(effects: &str) -> Vec<&'static str> {
    effects.split('&').map(|e| VARIANT_EFFECTS.iter().find(|(_, s)| *s == e).map(|(n, _)| *n).unwrap_or("OTHER")).collect()
}

fn coding_effect(s: &str) -> Result<String> {
    if CODING_EFFECTS.contains(&s) {
        Ok(s.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.CodingEffect.{s}")
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AllelicDepth {
    pub total_read_count: i32,
    pub allele_read_count: i32,
}

impl AllelicDepth {
    /// AllelicDepth.containsAllelicDepth
    pub fn contains_allelic_depth(g: Option<&Genotype>) -> bool {
        g.is_some_and(|g| g.ad.as_ref().is_some_and(|ad| ad.len() > 1))
    }
    /// AllelicDepth.fromGenotype (the genotype has AD)
    pub fn from_genotype(g: &Genotype) -> AllelicDepth {
        let ad = g.ad.as_ref().expect("genotype with AD");
        let sum = ad.iter().fold(0i32, |a, b| a.wrapping_add(*b));
        let total = match g.dp {
            Some(dp) => dp.max(sum),
            None => sum,
        };
        AllelicDepth { total_read_count: total, allele_read_count: ad[1] }
    }
}

/// VariantImpact (VariantImpactSerialiser.fromVariantContext)
#[derive(Clone, Debug, PartialEq)]
pub struct VariantImpact {
    pub gene_name: String,
    pub canonical_transcript: String,
    pub canonical_effect: String,
    pub canonical_coding_effect: String,
    pub canonical_hgvs_coding: String,
    pub canonical_hgvs_protein: String,
    pub canonical_splice_region: bool,
    pub other_reportable_effects: String,
    pub worst_coding_effect: String,
    pub genes_affected: i32,
}

impl VariantImpact {
    pub fn from_record(rec: &VcfRecord) -> Result<VariantImpact> {
        let v = rec.info.as_string_list("IMPACT", "");
        if v.len() != 10 {
            return Ok(VariantImpact {
                gene_name: String::new(),
                canonical_transcript: String::new(),
                canonical_effect: String::new(),
                canonical_coding_effect: "UNDEFINED".into(),
                canonical_hgvs_coding: String::new(),
                canonical_hgvs_protein: String::new(),
                canonical_splice_region: false,
                other_reportable_effects: String::new(),
                worst_coding_effect: "UNDEFINED".into(),
                genes_affected: 0,
            });
        }
        Ok(VariantImpact {
            gene_name: v[0].clone(),
            canonical_transcript: v[1].clone(),
            canonical_effect: v[2].clone(),
            canonical_coding_effect: coding_effect(&v[3])?,
            canonical_splice_region: crate::java_text::parse_bool(&v[4]),
            canonical_hgvs_coding: v[5].clone(),
            canonical_hgvs_protein: v[6].clone(),
            other_reportable_effects: v[7].clone(),
            worst_coding_effect: coding_effect(&v[8])?,
            genes_affected: crate::java_text::parse_int(&v[9])
                .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{}\"", v[9]))?,
        })
    }
}

/// VariantTranscriptImpact (PAVE_TI)
#[derive(Clone, Debug, PartialEq)]
pub struct VariantTranscriptImpact {
    pub gene_id: String,
    pub gene_name: String,
    pub transcript: String,
    pub effects: String,
    pub splice_region: bool,
    pub hgvs_coding: String,
    pub hgvs_protein: String,
    pub ref_seq_id: String,
    pub affected_exon: i32,
    pub affected_codon: i32,
}

impl VariantTranscriptImpact {
    /// fromVcfData: `|`-separated, the RefSeq id, exon and codon from v1.8 (10 items or more)
    pub fn from_vcf_data(data: &str) -> Result<VariantTranscriptImpact> {
        let items: Vec<&str> = data.split('|').collect();
        if items.len() < 7 {
            bail!("java.lang.ArrayIndexOutOfBoundsException: PAVE_TI item({data})");
        }
        let int = |s: &str| crate::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""));
        let (ref_seq_id, affected_exon, affected_codon) =
            if items.len() >= 10 { (items[7].to_string(), int(items[8])?, int(items[9])?) } else { (String::new(), 0, 0) };
        Ok(VariantTranscriptImpact {
            gene_id: items[0].into(),
            gene_name: items[1].into(),
            transcript: items[2].into(),
            effects: items[3].into(),
            splice_region: crate::java_text::parse_bool(items[4]),
            hgvs_coding: items[5].into(),
            hgvs_protein: items[6].into(),
            ref_seq_id,
            affected_exon,
            affected_codon,
        })
    }

    /// fromVariantContext
    pub fn from_record(rec: &VcfRecord) -> Result<Vec<VariantTranscriptImpact>> {
        if !rec.info.has("PAVE_TI") {
            return Ok(Vec::new());
        }
        rec.info.as_string_list("PAVE_TI", "").iter().map(|s| Self::from_vcf_data(s)).collect()
    }
}

/// AltTranscriptReportableInfo
#[derive(Clone, Debug, PartialEq)]
pub struct AltTranscriptReportableInfo {
    pub gene_name: String,
    pub trans_name: String,
    pub hgvs_coding: String,
    pub hgvs_protein: String,
    pub effects: String,
    pub effect: String,
}

impl AltTranscriptReportableInfo {
    /// parse: 6 items (with the gene) or 5 (without); None for fewer than 5
    pub fn parse(info: &str) -> Result<Option<AltTranscriptReportableInfo>> {
        let v: Vec<&str> = info.split('|').collect();
        if v.len() < 5 {
            return Ok(None);
        }
        let (gene_name, rest) = if v.len() == 6 { (v[0].to_string(), &v[1..]) } else { (String::new(), &v[..]) };
        Ok(Some(AltTranscriptReportableInfo {
            gene_name,
            trans_name: rest[0].into(),
            hgvs_coding: rest[1].into(),
            hgvs_protein: rest[2].into(),
            effects: rest[3].into(),
            effect: coding_effect(rest[4])?,
        }))
    }

    /// parseAltTranscriptInfo: the "--"-separated items
    pub fn parse_all(other_reportable_effects: &str) -> Result<Vec<AltTranscriptReportableInfo>> {
        let mut out = Vec::new();
        for s in other_reportable_effects.split("--") {
            if let Some(i) = Self::parse(s)? {
                out.push(i);
            }
        }
        Ok(out)
    }
}

/// Pathogenicity.fromClinvarAnnotation
pub fn pathogenicity_from_clinvar(sig: &str, conf: &str) -> &'static str {
    let p = |s: &str| s.contains("Pathogenic");
    let lp = |s: &str| s.contains("Likely_pathogenic");
    let b = |s: &str| s.contains("Benign");
    let lb = |s: &str| s.contains("Likely_benign");
    if sig.is_empty() {
        return "UNKNOWN";
    }
    if p(sig) {
        return "PATHOGENIC";
    }
    if lp(sig) {
        return "LIKELY_PATHOGENIC";
    }
    if b(sig) {
        return "BENIGN";
    }
    if lb(sig) {
        return "LIKELY_BENIGN";
    }
    if sig.contains("Conflicting") {
        if p(conf) && !b(conf) && !lb(conf) {
            return "PATHOGENIC";
        }
        if lp(conf) && !p(conf) && !b(conf) && !lb(conf) {
            return "LIKELY_PATHOGENIC";
        }
        if b(conf) && !p(conf) && !lp(conf) {
            return "BENIGN";
        }
        if lb(conf) && !p(conf) && !lp(conf) && !b(conf) {
            return "LIKELY_BENIGN";
        }
        if (p(conf) || lp(conf)) && (b(conf) || lb(conf)) {
            return "CONFLICTING";
        }
    }
    "UNKNOWN"
}

/// PathogenicSummaryFactory.fromContext: (clinvar info, Pathogenicity name)
pub fn pathogenic_summary(rec: &VcfRecord) -> (String, &'static str) {
    let sig = rec.info.as_string_list("CLNSIG", "").join(",");
    let conf = rec.info.as_string_list("CLNSIGCONF", "").join(",");
    let info = if !conf.is_empty() { conf.clone() } else { sig.clone() };
    let blacklisted = rec.info.as_bool("BLACKLIST_BED", false) || rec.info.as_bool("BLACKLIST_VCF", false);
    let status = if blacklisted { "BENIGN_BLACKLIST" } else { pathogenicity_from_clinvar(&sig, &conf) };
    (info, status)
}

/// Pathogenicity.isPathogenic
pub fn is_pathogenic(status: &str) -> bool {
    status == "PATHOGENIC" || status == "LIKELY_PATHOGENIC"
}

/// SmallVariant as SmallVariantFactory builds it
#[derive(Clone, Debug, PartialEq)]
pub struct SmallVariant {
    pub chromosome: String,
    pub position: i32,
    /// VariantType name: MNP, SNP, INDEL, UNDEFINED
    pub variant_type: &'static str,
    pub gene: String,
    pub ref_bases: String,
    pub alt: String,
    pub allelic_depth: AllelicDepth,
    pub reference_depth: Option<AllelicDepth>,
    pub rna_depth: Option<AllelicDepth>,
    pub canonical_transcript: String,
    pub canonical_effect: String,
    pub canonical_coding_effect: String,
    pub canonical_hgvs_coding_impact: String,
    pub canonical_hgvs_protein_impact: String,
    pub canonical_affected_codon: i32,
    pub canonical_affected_exon: i32,
    pub qual: f64,
    pub mappability: f64,
    pub filter: String,
    pub genes_affected: i32,
    pub splice_region: bool,
    pub other_reported_effects: String,
    pub worst_coding_effect: String,
    /// HotspotType name
    pub hotspot: &'static str,
    pub adjusted_copy_number: f64,
    pub adjusted_vaf: f64,
    pub minor_allele_copy_number: f64,
    pub variant_copy_number: f64,
    pub biallelic: bool,
    pub biallelic_probability: f64,
    pub reported: bool,
    /// GenotypeStatus name of the reference: HOM_REF, HET, HOM_ALT, UNKNOWN
    pub genotype_status: &'static str,
    pub germline_status: String,
    pub trinucleotide_context: String,
    pub microhomology: String,
    pub repeat_sequence: String,
    pub repeat_count: i32,
    pub tier: String,
    pub reportable_transcripts: Vec<String>,
    pub clinvar_info: String,
    pub subclonal_likelihood: f64,
    pub kataegis: String,
    pub gnomad_frequency: f64,
    /// SomaticLikelihood name
    pub somatic_likelihood: String,
    pub pathogenic: bool,
    pub pathogenicity: String,
    pub local_phase_sets: Option<Vec<i32>>,
}

/// GenotypeStatus.fromGenotype
fn genotype_status(g: Option<&Genotype>) -> &'static str {
    match g.map(|g| g.genotype_type()) {
        Some(GenotypeType::Het) => "HET",
        Some(GenotypeType::HomRef) => "HOM_REF",
        Some(GenotypeType::HomVar) => "HOM_ALT",
        _ => "UNKNOWN",
    }
}

const GERMLINE_STATUSES: [&str; 9] =
    ["HOM_DELETION", "HET_DELETION", "AMPLIFICATION", "NOISE", "DIPLOID", "LIKELY_DIPLOID", "EXCLUDED", "CENTROMETIC", "UNKNOWN"];

/// HumanChromosome.contains
pub fn is_human_chromosome(chromosome: &str) -> bool {
    let c = chromosome.strip_prefix("chr").unwrap_or(chromosome);
    if !c.is_empty() && c.bytes().all(|b| b.is_ascii_digit()) {
        return c.parse::<i64>().is_ok_and(|n| (1..=22).contains(&n));
    }
    c == "X" || c == "Y"
}

/// The filters SmallVariantFactory applies before building a variant
pub fn passes_factory_filters(rec: &VcfRecord, pass_only: bool) -> bool {
    if pass_only && rec.is_filtered() {
        return false;
    }
    if !is_human_chromosome(&rec.chrom) {
        return false;
    }
    match rec.info.get("NT") {
        None => true,
        Some(v) => v.java_string() == "ref",
    }
}

/// SmallVariantFactory.createVariant (the record passed the filters): None where the tumor has no
/// allelic depth or no reads
pub fn create_variant(
    rec: &mut VcfRecord,
    header: &VcfHeader,
    tumor: &str,
    reference: Option<&str>,
    rna: Option<&str>,
) -> Result<Option<SmallVariant>> {
    rec.decode(header)?;
    let tumor_g = rec.decoded_genotype(tumor);
    if !AllelicDepth::contains_allelic_depth(tumor_g) {
        return Ok(None);
    }
    let tumor_depth = AllelicDepth::from_genotype(tumor_g.unwrap());
    if tumor_depth.total_read_count <= 0 {
        return Ok(None);
    }
    let reference_depth = reference
        .and_then(|r| rec.decoded_genotype(r))
        .filter(|g| AllelicDepth::contains_allelic_depth(Some(g)))
        .map(AllelicDepth::from_genotype);
    let rna_depth = rna
        .and_then(|r| rec.decoded_genotype(r))
        .filter(|g| AllelicDepth::contains_allelic_depth(Some(g)))
        .map(AllelicDepth::from_genotype);

    let variant_type = match rec.variant_type() {
        VariantType::Mnp => "MNP",
        VariantType::Snp => "SNP",
        VariantType::Indel => "INDEL",
        _ => "UNDEFINED",
    };
    let impact = VariantImpact::from_record(rec)?;
    let genotype_status = match reference {
        Some(r) => genotype_status(rec.decoded_genotype(r)),
        None => "UNKNOWN",
    };
    let info = &rec.info;
    let pathogenic = if info.has("CLNSIG") { Some(pathogenic_summary(rec)) } else { None };

    let transcript_impacts = VariantTranscriptImpact::from_record(rec)?;
    let (mut canonical_affected_codon, mut canonical_affected_exon) = (0, 0);
    if let Some(t) = transcript_impacts.iter().find(|t| t.transcript == impact.canonical_transcript) {
        canonical_affected_exon = t.affected_exon;
        canonical_affected_codon = t.affected_codon;
    }

    let biallelic = info.as_bool("BIALLELIC", false);
    let hotspot = if info.as_bool("HOTSPOT", false) {
        "HOTSPOT"
    } else if info.as_bool("NEAR_HOTSPOT", false) {
        "NEAR_HOTSPOT"
    } else {
        "NON_HOTSPOT"
    };
    let reportable_transcripts = if info.has("REPORTABLE_TRANSCRIPTS") {
        info.as_string("REPORTABLE_TRANSCRIPTS", "").split('|').map(String::from).collect()
    } else {
        Vec::new()
    };
    let germline_status = {
        let s = info.as_string("PURPLE_GERMLINE", "UNKNOWN");
        if GERMLINE_STATUSES.contains(&s.as_str()) { s } else { "UNKNOWN".to_string() }
    };
    let somatic_likelihood = {
        let s = info.as_string("SOM_LH", "UNKNOWN");
        if !["UNKNOWN", "HIGH", "MEDIUM", "LOW"].contains(&s.as_str()) {
            bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.SomaticLikelihood.{s}");
        }
        s
    };
    let tier = match info.as_string("TIER", "UNKNOWN").as_str() {
        t @ ("HOTSPOT" | "PANEL" | "HIGH_CONFIDENCE" | "LOW_CONFIDENCE") => t.to_string(),
        _ => "UNKNOWN".to_string(),
    };

    Ok(Some(SmallVariant {
        chromosome: rec.chrom.clone(),
        position: rec.pos,
        variant_type,
        gene: impact.gene_name.clone(),
        ref_bases: rec.ref_allele.clone(),
        alt: rec.alts.join(","),
        allelic_depth: tumor_depth,
        reference_depth,
        rna_depth,
        canonical_transcript: impact.canonical_transcript.clone(),
        canonical_effect: impact.canonical_effect.clone(),
        canonical_coding_effect: impact.canonical_coding_effect.clone(),
        canonical_hgvs_coding_impact: impact.canonical_hgvs_coding.clone(),
        canonical_hgvs_protein_impact: impact.canonical_hgvs_protein.clone(),
        canonical_affected_codon,
        canonical_affected_exon,
        qual: rec.phred_scaled_qual(),
        mappability: info.as_double("MAPPABILITY", 0.0),
        filter: rec.display_filter(),
        genes_affected: impact.genes_affected,
        splice_region: impact.canonical_splice_region,
        other_reported_effects: impact.other_reportable_effects.clone(),
        worst_coding_effect: impact.worst_coding_effect.clone(),
        hotspot,
        adjusted_copy_number: info.as_double("PURPLE_CN", 0.0),
        adjusted_vaf: info.as_double("PURPLE_AF", 0.0),
        minor_allele_copy_number: info.as_double("PURPLE_MACN", 0.0),
        variant_copy_number: info.as_double("PURPLE_VCN", 0.0),
        biallelic,
        biallelic_probability: info.as_double("BIALLELIC_PROB", if biallelic { 1.0 } else { 0.0 }),
        reported: info.as_bool("REPORTED", false),
        genotype_status,
        germline_status,
        trinucleotide_context: info.as_string("TNC", ""),
        microhomology: info.as_string("MH", ""),
        repeat_sequence: info.as_string("REP_S", ""),
        repeat_count: info.as_int("REP_C", 0),
        tier,
        reportable_transcripts,
        clinvar_info: info.as_string("CLNSIG", ""),
        subclonal_likelihood: info.as_double("SUBCL", 0.0),
        kataegis: info.as_string("KT", ""),
        gnomad_frequency: info.as_double("GND_FREQ", 0.0),
        somatic_likelihood,
        pathogenic: pathogenic.as_ref().is_some_and(|p| is_pathogenic(p.1)),
        pathogenicity: pathogenic.map(|p| p.1.to_string()).unwrap_or_default(),
        local_phase_sets: if info.has("LPS") { Some(info.as_int_list("LPS", 0)) } else { None },
    }))
}

/// SmallVariantFactory.fromVCFFile(tumor, reference, rna, vcfFile) with `pass_only` for passOnlyInstance()
pub fn read_small_variants(
    path: &str,
    tumor: &str,
    reference: Option<&str>,
    rna: Option<&str>,
    pass_only: bool,
) -> Result<Vec<SmallVariant>> {
    let mut reader = VcfReader::open(path)?;
    let header = reader.header.clone();
    let has = |s: &str| header.samples.iter().any(|x| x == s);
    if !has(tumor) {
        bail!("java.lang.IllegalArgumentException: Sample {tumor} not found in vcf file {path}");
    }
    if let Some(r) = reference {
        if !has(r) {
            bail!("java.lang.IllegalArgumentException: Sample {r} not found in vcf file {path}");
        }
    }
    if let Some(r) = rna {
        if !has(r) {
            bail!("java.lang.IllegalArgumentException: Sample {r} not found in vcf file {path}");
        }
    }
    if !header.has_format("AD") {
        bail!("java.lang.IllegalArgumentException: Allelic depths is a required format field in vcf file {path}");
    }
    let mut out = Vec::new();
    while let Some(mut rec) = reader.next_record()? {
        if !passes_factory_filters(&rec, pass_only) {
            continue;
        }
        if let Some(v) = create_variant(&mut rec, &header, tumor, reference, rna)? {
            out.push(v);
        }
    }
    Ok(out)
}

/// The sample names of a VCF's header (VcfFileReader.genotypeOrdinals' keys)
pub fn vcf_samples(path: &str) -> Result<Vec<String>> {
    Ok(VcfReader::open(path)?.header.samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn effects_and_alt_transcripts() {
        assert_eq!(effects_to_list("missense_variant&splice_region_variant"), vec!["MISSENSE", "OTHER"]);
        assert_eq!(effects_to_list("stop_gained"), vec!["STOP_GAINED"]);
        let alts = AltTranscriptReportableInfo::parse_all("GENE|ENST2|c.1A>T|p.Met1?|start_lost|NONSENSE_OR_FRAMESHIFT--ENST3|c.2|p.X|x|MISSENSE").unwrap();
        assert_eq!(alts.len(), 2);
        assert_eq!(alts[0].gene_name, "GENE");
        assert_eq!(alts[1].gene_name, "");
        assert_eq!(alts[1].trans_name, "ENST3");
        assert!(AltTranscriptReportableInfo::parse_all("").unwrap().is_empty());
        assert!(AltTranscriptReportableInfo::parse("a|b|c|d|NOT_AN_EFFECT").is_err());
    }

    #[test]
    fn human_chromosomes() {
        assert!(is_human_chromosome("chr1") && is_human_chromosome("22") && is_human_chromosome("chrX"));
        assert!(!is_human_chromosome("chrM") && !is_human_chromosome("23") && !is_human_chromosome("chr1_alt") && !is_human_chromosome("0"));
    }

    #[test]
    fn transcript_impacts() {
        let t = VariantTranscriptImpact::from_vcf_data("ENSG1|G|ENST1|missense_variant|false|c.1A>T|p.K1N|NM_1|3|45").unwrap();
        assert_eq!((t.affected_exon, t.affected_codon), (3, 45));
        let t = VariantTranscriptImpact::from_vcf_data("ENSG1|G|ENST1|missense_variant|true|c.1A>T|p.K1N").unwrap();
        assert!(t.splice_region && t.affected_exon == 0);
    }
}

/// CodingEffect.effect(VariantEffect), by enum names
pub fn coding_effect_of(variant_effect: &str) -> &'static str {
    match variant_effect {
        "FRAMESHIFT" | "STOP_GAINED" | "STOP_LOST" | "START_LOST" => "NONSENSE_OR_FRAMESHIFT",
        "SPLICE_ACCEPTOR" | "SPLICE_DONOR" => "SPLICE",
        "MISSENSE" | "INFRAME_DELETION" | "INFRAME_INSERTION" | "PHASED_INFRAME_DELETION" | "PHASED_INFRAME_INSERTION" | "PHASED_MISSENSE" => {
            "MISSENSE"
        }
        "SYNONYMOUS" | "PHASED_SYNONYMOUS" => "SYNONYMOUS",
        _ => "NONE",
    }
}

/// ORANGE's PurpleConversion.determineCodingEffect: the most severe coding effect of the effects
pub fn determine_coding_effect(variant_effects: &[&str]) -> &'static str {
    let simplified: Vec<&str> = variant_effects.iter().map(|e| coding_effect_of(e)).collect();
    for e in ["NONSENSE_OR_FRAMESHIFT", "SPLICE", "MISSENSE", "SYNONYMOUS"] {
        if simplified.contains(&e) {
            return e;
        }
    }
    "NONE"
}
