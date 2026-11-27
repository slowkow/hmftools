//! The ORANGE data model (`orange-datamodel`, com.hartwig.hmftools.datamodel) as the JSON report holds it.
//! Each type writes its attributes in the order of its interface, as the Immutables Gson adapters do; a
//! @Nullable attribute that is absent is written as null. Enum values are kept as their names; a
//! Set<Enum> is an EnumSet (ordinal order, no repeats): `enum_set` builds it from the enum's names.

use crate::json::{JsonWriter, ToJson};

/// The elements of an EnumSet built from `items`: ordered as `order`, without repeats
pub fn enum_set<S: AsRef<str>>(items: &[S], order: &[&str]) -> Vec<String> {
    order.iter().filter(|o| items.iter().any(|i| i.as_ref() == **o)).map(|o| o.to_string()).collect()
}

// datamodel enums whose order matters (EnumSet attributes)
pub const PURPLE_QC_STATUSES: [&str; 9] = [
    "PASS",
    "WARN_DELETED_GENES",
    "WARN_HIGH_COPY_NUMBER_NOISE",
    "WARN_GENDER_MISMATCH",
    "WARN_LOW_PURITY",
    "WARN_TINC",
    "FAIL_CONTAMINATION",
    "FAIL_NO_TUMOR",
    "FAIL_TINC",
];
pub const PURPLE_GERMLINE_ABERRATIONS: [&str; 9] =
    ["NONE", "MOSAIC_X", "KLINEFELTER", "XYY", "TRISOMY_X", "TRISOMY_13", "TRISOMY_15", "TRISOMY_18", "TRISOMY_21"];
pub const PURPLE_VARIANT_EFFECTS: [&str; 20] = [
    "STOP_GAINED",
    "STOP_LOST",
    "START_LOST",
    "FRAMESHIFT",
    "SPLICE_ACCEPTOR",
    "SPLICE_DONOR",
    "INFRAME_INSERTION",
    "INFRAME_DELETION",
    "MISSENSE",
    "PHASED_MISSENSE",
    "PHASED_INFRAME_INSERTION",
    "PHASED_INFRAME_DELETION",
    "SYNONYMOUS",
    "PHASED_SYNONYMOUS",
    "INTRONIC",
    "FIVE_PRIME_UTR",
    "THREE_PRIME_UTR",
    "UPSTREAM_GENE",
    "NON_CODING_TRANSCRIPT",
    "OTHER",
];
pub const RNA_QC_STATUSES: [&str; 4] = ["PASS", "FAIL_LOW_COVERAGE", "WARN_DUPLICATE_RATE", "WARN_SPLICED_GENE_COVERAGE"];
/// PurpleGermlineStatus (fewer values than hmf-common's GermlineStatus)
pub const PURPLE_GERMLINE_STATUSES: [&str; 6] = ["HOM_DELETION", "HET_DELETION", "AMPLIFICATION", "NOISE", "DIPLOID", "UNKNOWN"];

// ---------------------------------------------------------------------------------------------
// orange

#[derive(Clone, Debug, PartialEq)]
pub struct OrangeDoidNode {
    pub doid: String,
    pub doid_term: Option<String>,
}

impl ToJson for OrangeDoidNode {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object().f_str("doid", &self.doid).f_opt_str("doidTerm", self.doid_term.as_deref()).end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct OrangePlot {
    pub filename: String,
    pub title: String,
}

impl ToJson for OrangePlot {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object().f_str("filename", &self.filename).f_str("title", &self.title).end_object();
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrangePlots {
    pub purple_input_circos_plot: String,
    pub purple_final_circos_plot: String,
    pub purple_clonality_plot: String,
    pub purple_copy_number_plot: String,
    pub purple_minor_allele_map_plot: String,
    pub purple_variant_copy_number_plot: String,
    pub purple_purity_range_plot: String,
    pub purple_rainfall_plot: String,
    pub linx_driver_plots: Vec<OrangePlot>,
    pub cuppa_summary_plot: Option<String>,
    pub q_see_plot: Option<String>,
}

impl ToJson for OrangePlots {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("purpleInputCircosPlot", &self.purple_input_circos_plot)
            .f_str("purpleFinalCircosPlot", &self.purple_final_circos_plot)
            .f_str("purpleClonalityPlot", &self.purple_clonality_plot)
            .f_str("purpleCopyNumberPlot", &self.purple_copy_number_plot)
            .f_str("purpleMinorAlleleMapPlot", &self.purple_minor_allele_map_plot)
            .f_str("purpleVariantCopyNumberPlot", &self.purple_variant_copy_number_plot)
            .f_str("purplePurityRangePlot", &self.purple_purity_range_plot)
            .f_str("purpleRainfallPlot", &self.purple_rainfall_plot)
            .f_list("linxDriverPlots", &self.linx_driver_plots)
            .f_opt_str("cuppaSummaryPlot", self.cuppa_summary_plot.as_deref())
            .f_opt_str("qSeePlot", self.q_see_plot.as_deref())
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct OrangeRecord {
    pub sample_id: String,
    pub reference_id: Option<String>,
    pub rna_id: Option<String>,
    /// LocalDate.toString (yyyy-MM-dd)
    pub sampling_date: String,
    /// ExperimentType name
    pub experiment_type: String,
    pub configured_primary_tumor: Vec<OrangeDoidNode>,
    pub ref_genome_version: String,
    pub pipeline_version: Option<String>,
    pub panel_name: Option<String>,
    pub sequencing_type: String,
    pub purple: PurpleRecord,
    pub linx: LinxRecord,
    pub isofox: Option<IsofoxRecord>,
    pub lilac: Option<LilacRecord>,
    pub immune_escape: ImmuneEscapeRecord,
    pub virus_interpreter: Option<VirusInterpreterData>,
    pub chord: Option<ChordRecord>,
    pub cuppa: Option<CuppaData>,
    pub peach: Option<Vec<PeachGenotype>>,
    pub sig_allocations: Option<Vec<SignatureAllocation>>,
    pub plots: OrangePlots,
}

impl ToJson for OrangeRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("sampleId", &self.sample_id)
            .f_opt_str("referenceId", self.reference_id.as_deref())
            .f_opt_str("rnaId", self.rna_id.as_deref())
            .f_str("samplingDate", &self.sampling_date)
            .f_str("experimentType", &self.experiment_type)
            .f_list("configuredPrimaryTumor", &self.configured_primary_tumor)
            .f_str("refGenomeVersion", &self.ref_genome_version)
            .f_opt_str("pipelineVersion", self.pipeline_version.as_deref())
            .f_opt_str("panelName", self.panel_name.as_deref())
            .f_str("sequencingType", &self.sequencing_type)
            .f_value("purple", &self.purple)
            .f_value("linx", &self.linx)
            .f_opt("isofox", self.isofox.as_ref())
            .f_opt("lilac", self.lilac.as_ref())
            .f_value("immuneEscape", &self.immune_escape)
            .f_opt("virusInterpreter", self.virus_interpreter.as_ref())
            .f_opt("chord", self.chord.as_ref())
            .f_opt("cuppa", self.cuppa.as_ref())
            .f_opt_list("peach", self.peach.as_deref())
            .f_opt_list("sigAllocations", self.sig_allocations.as_deref())
            .f_value("plots", &self.plots)
            .end_object();
    }
}

// ---------------------------------------------------------------------------------------------
// purple

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleQC {
    pub status: Vec<String>,
    pub germline_aberrations: Vec<String>,
    pub amber_mean_depth: i32,
    pub contamination: f64,
    pub total_copy_number_segments: i32,
    pub unsupported_copy_number_segments: i32,
    pub deleted_genes: i32,
}

impl ToJson for PurpleQC {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str_list("status", &self.status)
            .f_str_list("germlineAberrations", &self.germline_aberrations)
            .f_int("amberMeanDepth", self.amber_mean_depth as i64)
            .f_double("contamination", self.contamination)
            .f_int("totalCopyNumberSegments", self.total_copy_number_segments as i64)
            .f_int("unsupportedCopyNumberSegments", self.unsupported_copy_number_segments as i64)
            .f_int("deletedGenes", self.deleted_genes as i64)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleFit {
    pub qc: PurpleQC,
    pub fitted_purity_method: String,
    pub purity: f64,
    pub min_purity: f64,
    pub max_purity: f64,
    pub ploidy: f64,
    pub min_ploidy: f64,
    pub max_ploidy: f64,
}

impl ToJson for PurpleFit {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_value("qc", &self.qc)
            .f_str("fittedPurityMethod", &self.fitted_purity_method)
            .f_double("purity", self.purity)
            .f_double("minPurity", self.min_purity)
            .f_double("maxPurity", self.max_purity)
            .f_double("ploidy", self.ploidy)
            .f_double("minPloidy", self.min_ploidy)
            .f_double("maxPloidy", self.max_ploidy)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleCharacteristics {
    pub whole_genome_duplication: bool,
    pub microsatellite_indels_per_mb: f64,
    pub microsatellite_status: String,
    pub tumor_mutational_burden_per_mb: f64,
    pub tumor_mutational_burden_status: String,
    pub tumor_mutational_load: i32,
    pub tumor_mutational_load_status: String,
    pub sv_tumor_mutational_burden: i32,
    pub loh_percentage: f64,
}

impl ToJson for PurpleCharacteristics {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_bool("wholeGenomeDuplication", self.whole_genome_duplication)
            .f_double("microsatelliteIndelsPerMb", self.microsatellite_indels_per_mb)
            .f_str("microsatelliteStatus", &self.microsatellite_status)
            .f_double("tumorMutationalBurdenPerMb", self.tumor_mutational_burden_per_mb)
            .f_str("tumorMutationalBurdenStatus", &self.tumor_mutational_burden_status)
            .f_int("tumorMutationalLoad", self.tumor_mutational_load as i64)
            .f_str("tumorMutationalLoadStatus", &self.tumor_mutational_load_status)
            .f_int("svTumorMutationalBurden", self.sv_tumor_mutational_burden as i64)
            .f_double("lohPercentage", self.loh_percentage)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleDriver {
    pub gene: String,
    pub transcript: String,
    /// PurpleDriverType name
    pub driver_type: String,
    pub driver_likelihood: f64,
    pub likelihood_method: String,
    pub reported_status: String,
    pub driver_interpretation: String,
    pub is_canonical: bool,
    pub category: String,
}

impl ToJson for PurpleDriver {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("gene", &self.gene)
            .f_str("transcript", &self.transcript)
            .f_str("type", &self.driver_type)
            .f_double("driverLikelihood", self.driver_likelihood)
            .f_str("likelihoodMethod", &self.likelihood_method)
            .f_str("reportedStatus", &self.reported_status)
            .f_str("driverInterpretation", &self.driver_interpretation)
            .f_bool("isCanonical", self.is_canonical)
            .f_str("category", &self.category)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AllelicDepth {
    pub total_read_count: i32,
    pub allele_read_count: i32,
}

impl ToJson for AllelicDepth {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_int("totalReadCount", self.total_read_count as i64)
            .f_int("alleleReadCount", self.allele_read_count as i64)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleTranscriptImpact {
    pub transcript: String,
    pub hgvs_coding_impact: String,
    pub hgvs_protein_impact: String,
    pub affected_codon: Option<i32>,
    pub affected_exon: Option<i32>,
    pub in_splice_region: bool,
    /// an EnumSet of PurpleVariantEffect names
    pub effects: Vec<String>,
    pub coding_effect: String,
    pub reported: bool,
}

impl ToJson for PurpleTranscriptImpact {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("transcript", &self.transcript)
            .f_str("hgvsCodingImpact", &self.hgvs_coding_impact)
            .f_str("hgvsProteinImpact", &self.hgvs_protein_impact)
            .f_opt_int("affectedCodon", self.affected_codon.map(i64::from))
            .f_opt_int("affectedExon", self.affected_exon.map(i64::from))
            .f_bool("inSpliceRegion", self.in_splice_region)
            .f_str_list("effects", &self.effects)
            .f_str("codingEffect", &self.coding_effect)
            .f_bool("reported", self.reported)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleVariant {
    pub variant_type: String,
    pub gene: String,
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt: String,
    pub driver_likelihood: f64,
    pub worst_coding_effect: String,
    pub canonical_impact: PurpleTranscriptImpact,
    pub other_impacts: Vec<PurpleTranscriptImpact>,
    pub hotspot: String,
    pub tumor_depth: AllelicDepth,
    pub rna_depth: Option<AllelicDepth>,
    pub adjusted_copy_number: f64,
    pub adjusted_vaf: f64,
    pub minor_allele_copy_number: f64,
    pub variant_copy_number: f64,
    pub biallelic: bool,
    pub biallelic_probability: Option<f64>,
    pub genotype_status: String,
    pub repeat_count: i32,
    pub subclonal_likelihood: f64,
    pub somatic_likelihood: Option<String>,
    pub local_phase_sets: Option<Vec<i32>>,
    pub clinvar_pathogenicity: Option<String>,
    pub gnomad_frequency: Option<f64>,
    pub plot_filename: Option<String>,
}

impl PurpleVariant {
    /// the @Value.Derived reported(): the canonical impact's, else any other impact's
    pub fn reported(&self) -> bool {
        self.canonical_impact.reported || self.other_impacts.iter().any(|i| i.reported)
    }
}

impl ToJson for PurpleVariant {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("type", &self.variant_type)
            .f_str("gene", &self.gene)
            .f_str("chromosome", &self.chromosome)
            .f_int("position", self.position as i64)
            .f_str("ref", &self.ref_bases)
            .f_str("alt", &self.alt)
            .f_double("driverLikelihood", self.driver_likelihood)
            .f_str("worstCodingEffect", &self.worst_coding_effect)
            .f_value("canonicalImpact", &self.canonical_impact)
            .f_list("otherImpacts", &self.other_impacts)
            .f_str("hotspot", &self.hotspot)
            .f_value("tumorDepth", &self.tumor_depth)
            .f_opt("rnaDepth", self.rna_depth.as_ref())
            .f_double("adjustedCopyNumber", self.adjusted_copy_number)
            .f_double("adjustedVAF", self.adjusted_vaf)
            .f_double("minorAlleleCopyNumber", self.minor_allele_copy_number)
            .f_double("variantCopyNumber", self.variant_copy_number)
            .f_bool("biallelic", self.biallelic)
            .f_opt_double("biallelicProbability", self.biallelic_probability)
            .f_str("genotypeStatus", &self.genotype_status)
            .f_int("repeatCount", self.repeat_count as i64)
            .f_double("subclonalLikelihood", self.subclonal_likelihood)
            .f_opt_str("somaticLikelihood", self.somatic_likelihood.as_deref());
        w.name("localPhaseSets");
        match &self.local_phase_sets {
            Some(v) => {
                w.begin_array();
                for i in v {
                    w.int(*i as i64);
                }
                w.end_array();
            }
            None => {
                w.null();
            }
        }
        w.f_opt_str("clinvarPathogenicity", self.clinvar_pathogenicity.as_deref())
            .f_opt_double("gnomadFrequency", self.gnomad_frequency)
            .f_opt_str("plotFilename", self.plot_filename.as_deref())
            .f_bool("reported", self.reported())
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GermlineAmpDelFields {
    pub germline_status: String,
    pub somatic_status: String,
    pub germline_min_copy_number: f64,
}

impl ToJson for GermlineAmpDelFields {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("germlineStatus", &self.germline_status)
            .f_str("somaticStatus", &self.somatic_status)
            .f_double("germlineMinCopyNumber", self.germline_min_copy_number)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleGainDeletion {
    pub driver: PurpleDriver,
    pub germline_amp_del_fields: Option<GermlineAmpDelFields>,
    pub chromosome: String,
    pub chromosome_band: String,
    pub gene_range: String,
    pub exon_start: Option<i32>,
    pub exon_end: Option<i32>,
    pub min_copy_number: f64,
    pub max_copy_number: f64,
    pub relative_copy_number: f64,
    pub arm_copy_number: f64,
    pub min_minor_allele_copies: f64,
    pub tpm: Option<f64>,
    pub tpm_percentile: Option<f64>,
    pub tpm_fold_change: Option<f64>,
    pub plot_filename: Option<String>,
}

impl ToJson for PurpleGainDeletion {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_value("driver", &self.driver)
            .f_opt("germlineAmpDelFields", self.germline_amp_del_fields.as_ref())
            .f_str("chromosome", &self.chromosome)
            .f_str("chromosomeBand", &self.chromosome_band)
            .f_str("geneRange", &self.gene_range)
            .f_opt_int("exonStart", self.exon_start.map(i64::from))
            .f_opt_int("exonEnd", self.exon_end.map(i64::from))
            .f_double("minCopyNumber", self.min_copy_number)
            .f_double("maxCopyNumber", self.max_copy_number)
            .f_double("relativeCopyNumber", self.relative_copy_number)
            .f_double("armCopyNumber", self.arm_copy_number)
            .f_double("minMinorAlleleCopies", self.min_minor_allele_copies)
            .f_opt_double("tpm", self.tpm)
            .f_opt_double("tpmPercentile", self.tpm_percentile)
            .f_opt_double("tpmFoldChange", self.tpm_fold_change)
            .f_opt_str("plotFilename", self.plot_filename.as_deref())
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleChrArmCopyNumber {
    pub chromosome: String,
    pub arm: String,
    pub cn_type: String,
    pub copy_number: f64,
    pub relative_copy_number: f64,
    pub driver_interpretation: String,
}

impl ToJson for PurpleChrArmCopyNumber {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("chromosome", &self.chromosome)
            .f_str("arm", &self.arm)
            .f_str("type", &self.cn_type)
            .f_double("copyNumber", self.copy_number)
            .f_double("relativeCopyNumber", self.relative_copy_number)
            .f_str("driverInterpretation", &self.driver_interpretation)
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct PurpleRecord {
    pub fit: PurpleFit,
    pub characteristics: PurpleCharacteristics,
    pub somatic_drivers: Vec<PurpleDriver>,
    pub germline_drivers: Option<Vec<PurpleDriver>>,
    pub somatic_variants: Vec<PurpleVariant>,
    pub germline_variants: Option<Vec<PurpleVariant>>,
    pub somatic_gains_dels: Vec<PurpleGainDeletion>,
    pub germline_gains_dels: Option<Vec<PurpleGainDeletion>>,
    pub arm_copy_number_abberations: Vec<PurpleChrArmCopyNumber>,
}

impl ToJson for PurpleRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_value("fit", &self.fit)
            .f_value("characteristics", &self.characteristics)
            .f_list("somaticDrivers", &self.somatic_drivers)
            .f_opt_list("germlineDrivers", self.germline_drivers.as_deref())
            .f_list("somaticVariants", &self.somatic_variants)
            .f_opt_list("germlineVariants", self.germline_variants.as_deref())
            .f_list("somaticGainsDels", &self.somatic_gains_dels)
            .f_opt_list("germlineGainsDels", self.germline_gains_dels.as_deref())
            .f_list("armCopyNumberAbberations", &self.arm_copy_number_abberations)
            .end_object();
    }
}

// ---------------------------------------------------------------------------------------------
// linx

#[derive(Clone, Debug, PartialEq)]
pub struct LinxFusion {
    pub gene_up: String,
    pub cluster_id: i32,
    pub context_up: String,
    pub transcript_up: String,
    pub gene_down: String,
    pub context_down: String,
    pub transcript_down: String,
    pub reported_type: String,
    pub phased: String,
    pub driver_interpretation: String,
    pub fused_exon_up: i32,
    pub fused_exon_down: i32,
    pub chain_links: i32,
    pub chain_terminated: bool,
    pub domains_kept: String,
    pub domains_lost: String,
    pub junction_copy_number: f64,
    pub rna_support: Option<AllelicDepth>,
    pub plot_filename: Option<String>,
}

impl ToJson for LinxFusion {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("geneUp", &self.gene_up)
            .f_int("clusterId", self.cluster_id as i64)
            .f_str("contextUp", &self.context_up)
            .f_str("transcriptUp", &self.transcript_up)
            .f_str("geneDown", &self.gene_down)
            .f_str("contextDown", &self.context_down)
            .f_str("transcriptDown", &self.transcript_down)
            .f_str("reportedType", &self.reported_type)
            .f_str("phased", &self.phased)
            .f_str("driverInterpretation", &self.driver_interpretation)
            .f_int("fusedExonUp", self.fused_exon_up as i64)
            .f_int("fusedExonDown", self.fused_exon_down as i64)
            .f_int("chainLinks", self.chain_links as i64)
            .f_bool("chainTerminated", self.chain_terminated)
            .f_str("domainsKept", &self.domains_kept)
            .f_str("domainsLost", &self.domains_lost)
            .f_double("junctionCopyNumber", self.junction_copy_number)
            .f_opt("rnaSupport", self.rna_support.as_ref())
            .f_opt_str("plotFilename", self.plot_filename.as_deref())
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LinxBreakend {
    pub id: i32,
    pub sv_id: i32,
    pub cluster_id: i32,
    pub gene: String,
    pub position: i32,
    pub orientation: i32,
    pub chromosome: String,
    /// @NotNull in the data model but null where the chromosome has no band (Immutables would throw)
    pub chromosome_band: String,
    pub transcript: String,
    pub is_canonical: bool,
    pub gene_orientation: String,
    pub disruptive: bool,
    pub reported_status: String,
    pub undisrupted_copy_number: f64,
    pub breakend_type: String,
    pub region_type: String,
    pub coding_type: String,
    pub next_splice_exon_rank: i32,
    pub exon_up: i32,
    pub exon_down: i32,
    pub junction_copy_number: f64,
    pub driver_type: String,
    pub driver_likelihood: f64,
    pub plot_filename: Option<String>,
}

impl ToJson for LinxBreakend {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_int("id", self.id as i64)
            .f_int("svId", self.sv_id as i64)
            .f_int("clusterId", self.cluster_id as i64)
            .f_str("gene", &self.gene)
            .f_int("position", self.position as i64)
            .f_int("orientation", self.orientation as i64)
            .f_str("chromosome", &self.chromosome)
            .f_str("chromosomeBand", &self.chromosome_band)
            .f_str("transcript", &self.transcript)
            .f_bool("isCanonical", self.is_canonical)
            .f_str("geneOrientation", &self.gene_orientation)
            .f_bool("disruptive", self.disruptive)
            .f_str("reportedStatus", &self.reported_status)
            .f_double("undisruptedCopyNumber", self.undisrupted_copy_number)
            .f_str("type", &self.breakend_type)
            .f_str("regionType", &self.region_type)
            .f_str("codingType", &self.coding_type)
            .f_int("nextSpliceExonRank", self.next_splice_exon_rank as i64)
            .f_int("exonUp", self.exon_up as i64)
            .f_int("exonDown", self.exon_down as i64)
            .f_double("junctionCopyNumber", self.junction_copy_number)
            .f_str("driverType", &self.driver_type)
            .f_double("driverLikelihood", self.driver_likelihood)
            .f_opt_str("plotFilename", self.plot_filename.as_deref())
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct LinxRecord {
    pub fusions: Vec<LinxFusion>,
    pub somatic_breakends: Vec<LinxBreakend>,
    pub germline_breakends: Option<Vec<LinxBreakend>>,
}

impl ToJson for LinxRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_list("fusions", &self.fusions)
            .f_list("somaticBreakends", &self.somatic_breakends)
            .f_opt_list("germlineBreakends", self.germline_breakends.as_deref())
            .end_object();
    }
}

// ---------------------------------------------------------------------------------------------
// isofox

#[derive(Clone, Debug, PartialEq)]
pub struct RnaStatistics {
    pub qc_status: Vec<String>,
    pub total_fragments: i64,
    pub duplicate_fragments: i64,
    pub spliced_fragment_perc: f64,
    pub unspliced_fragment_perc: f64,
    pub alt_fragment_perc: f64,
    pub chimeric_fragment_perc: f64,
}

impl ToJson for RnaStatistics {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str_list("qcStatus", &self.qc_status)
            .f_int("totalFragments", self.total_fragments)
            .f_int("duplicateFragments", self.duplicate_fragments)
            .f_double("splicedFragmentPerc", self.spliced_fragment_perc)
            .f_double("unsplicedFragmentPerc", self.unspliced_fragment_perc)
            .f_double("altFragmentPerc", self.alt_fragment_perc)
            .f_double("chimericFragmentPerc", self.chimeric_fragment_perc)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GeneExpression {
    pub gene: String,
    pub tpm: f64,
    pub median_tpm_cohort: f64,
    pub percentile_cohort: f64,
    pub median_tpm_cancer: Option<f64>,
    pub percentile_cancer: Option<f64>,
}

impl ToJson for GeneExpression {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("gene", &self.gene)
            .f_double("tpm", self.tpm)
            .f_double("medianTpmCohort", self.median_tpm_cohort)
            .f_double("percentileCohort", self.percentile_cohort)
            .f_opt_double("medianTpmCancer", self.median_tpm_cancer)
            .f_opt_double("percentileCancer", self.percentile_cancer)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RnaFusion {
    pub gene_start: Option<String>,
    pub gene_end: Option<String>,
    pub chromosome_start: String,
    pub chromosome_end: String,
    pub position_start: i32,
    pub position_end: i32,
    pub junction_type_start: String,
    pub junction_type_end: String,
    pub known_type: String,
    pub sv_type: String,
    pub split_fragments: i32,
    pub realigned_frags: i32,
    pub discordant_frags: i32,
    pub depth_start: i32,
    pub depth_end: i32,
    pub cohort_frequency: i32,
}

impl ToJson for RnaFusion {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_opt_str("geneStart", self.gene_start.as_deref())
            .f_opt_str("geneEnd", self.gene_end.as_deref())
            .f_str("chromosomeStart", &self.chromosome_start)
            .f_str("chromosomeEnd", &self.chromosome_end)
            .f_int("positionStart", self.position_start as i64)
            .f_int("positionEnd", self.position_end as i64)
            .f_str("junctionTypeStart", &self.junction_type_start)
            .f_str("junctionTypeEnd", &self.junction_type_end)
            .f_str("knownType", &self.known_type)
            .f_str("svType", &self.sv_type)
            .f_int("splitFragments", self.split_fragments as i64)
            .f_int("realignedFrags", self.realigned_frags as i64)
            .f_int("discordantFrags", self.discordant_frags as i64)
            .f_int("depthStart", self.depth_start as i64)
            .f_int("depthEnd", self.depth_end as i64)
            .f_int("cohortFrequency", self.cohort_frequency as i64)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct NovelSpliceJunction {
    pub gene: String,
    pub chromosome: String,
    pub junction_start: i32,
    pub junction_end: i32,
    pub junction_type: String,
    pub exon_start: i32,
    pub exon_end: i32,
    pub fragment_count: i32,
    pub depth_start: i32,
    pub depth_end: i32,
    pub region_start: String,
    pub region_end: String,
    pub cohort_frequency: i32,
}

impl ToJson for NovelSpliceJunction {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("gene", &self.gene)
            .f_str("chromosome", &self.chromosome)
            .f_int("junctionStart", self.junction_start as i64)
            .f_int("junctionEnd", self.junction_end as i64)
            .f_str("type", &self.junction_type)
            .f_int("exonStart", self.exon_start as i64)
            .f_int("exonEnd", self.exon_end as i64)
            .f_int("fragmentCount", self.fragment_count as i64)
            .f_int("depthStart", self.depth_start as i64)
            .f_int("depthEnd", self.depth_end as i64)
            .f_str("regionStart", &self.region_start)
            .f_str("regionEnd", &self.region_end)
            .f_int("cohortFrequency", self.cohort_frequency as i64)
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct IsofoxRecord {
    pub summary: RnaStatistics,
    pub high_expression_genes: Vec<GeneExpression>,
    pub low_expression_genes: Vec<GeneExpression>,
    pub fusions: Vec<RnaFusion>,
    pub novel_splice_junctions: Vec<NovelSpliceJunction>,
}

impl ToJson for IsofoxRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_value("summary", &self.summary)
            .f_list("highExpressionGenes", &self.high_expression_genes)
            .f_list("lowExpressionGenes", &self.low_expression_genes)
            .f_list("fusions", &self.fusions)
            .f_list("novelSpliceJunctions", &self.novel_splice_junctions)
            .end_object();
    }
}

// ---------------------------------------------------------------------------------------------
// lilac, immune escape, virus, chord, cuppa, peach, sigs

#[derive(Clone, Debug, PartialEq)]
pub struct LilacAllele {
    pub gene_class: String,
    pub allele: String,
    pub qc_status: String,
    pub tumor_copy_number: f64,
    pub ref_fragments: Option<i32>,
    pub tumor_fragments: i32,
    pub rna_fragments: Option<i32>,
    pub somatic_missense: f64,
    pub somatic_nonsense_or_frameshift: f64,
    pub somatic_splice: f64,
    pub somatic_synonymous: f64,
    pub somatic_inframe_indel: f64,
}

impl ToJson for LilacAllele {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("geneClass", &self.gene_class)
            .f_str("allele", &self.allele)
            .f_str("qcStatus", &self.qc_status)
            .f_double("tumorCopyNumber", self.tumor_copy_number)
            .f_opt_int("refFragments", self.ref_fragments.map(i64::from))
            .f_int("tumorFragments", self.tumor_fragments as i64)
            .f_opt_int("rnaFragments", self.rna_fragments.map(i64::from))
            .f_double("somaticMissense", self.somatic_missense)
            .f_double("somaticNonsenseOrFrameshift", self.somatic_nonsense_or_frameshift)
            .f_double("somaticSplice", self.somatic_splice)
            .f_double("somaticSynonymous", self.somatic_synonymous)
            .f_double("somaticInframeIndel", self.somatic_inframe_indel)
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct LilacRecord {
    pub alleles: Vec<LilacAllele>,
}

impl ToJson for LilacRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object().f_list("alleles", &self.alleles).end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImmuneEscapeRecord {
    pub has_hla_escape: bool,
    pub has_antigen_presentation_pathway_escape: bool,
    pub has_ifn_gamma_pathway_escape: bool,
    pub has_pdl1_overexpression_escape: bool,
    pub has_cd58_inactivation_escape: bool,
    pub has_epigenetic_setdb1_escape: bool,
}

impl ToJson for ImmuneEscapeRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_bool("hasHlaEscape", self.has_hla_escape)
            .f_bool("hasAntigenPresentationPathwayEscape", self.has_antigen_presentation_pathway_escape)
            .f_bool("hasIFNGammaPathwayEscape", self.has_ifn_gamma_pathway_escape)
            .f_bool("hasPDL1OverexpressionEscape", self.has_pdl1_overexpression_escape)
            .f_bool("hasCD58InactivationEscape", self.has_cd58_inactivation_escape)
            .f_bool("hasEpigeneticSETDB1Escape", self.has_epigenetic_setdb1_escape)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VirusInterpreterEntry {
    pub name: String,
    pub qc_status: String,
    pub integrations: i32,
    pub interpretation: Option<String>,
    pub percentage_covered: f64,
    pub mean_coverage: f64,
    pub expected_clonal_coverage: Option<f64>,
    pub reported: bool,
    pub driver_interpretation: String,
}

impl ToJson for VirusInterpreterEntry {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("name", &self.name)
            .f_str("qcStatus", &self.qc_status)
            .f_int("integrations", self.integrations as i64)
            .f_opt_str("interpretation", self.interpretation.as_deref())
            .f_double("percentageCovered", self.percentage_covered)
            .f_double("meanCoverage", self.mean_coverage)
            .f_opt_double("expectedClonalCoverage", self.expected_clonal_coverage)
            .f_bool("reported", self.reported)
            .f_str("driverInterpretation", &self.driver_interpretation)
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct VirusInterpreterData {
    pub all_viruses: Vec<VirusInterpreterEntry>,
    pub reportable_viruses: Vec<VirusInterpreterEntry>,
}

impl ToJson for VirusInterpreterData {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object().f_list("allViruses", &self.all_viruses).f_list("reportableViruses", &self.reportable_viruses).end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChordRecord {
    pub brca1_value: f64,
    pub brca2_value: f64,
    pub hrd_value: f64,
    pub hr_status: String,
    pub hrd_type: String,
}

impl ToJson for ChordRecord {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_double("BRCA1Value", self.brca1_value)
            .f_double("BRCA2Value", self.brca2_value)
            .f_double("hrdValue", self.hrd_value)
            .f_str("hrStatus", &self.hr_status)
            .f_str("hrdType", &self.hrd_type)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CuppaPrediction {
    pub cancer_type: String,
    pub likelihood: f64,
    pub snv_pairwise_classifier: Option<f64>,
    pub genomic_position_classifier: Option<f64>,
    pub feature_classifier: Option<f64>,
    pub alt_sj_cohort_classifier: Option<f64>,
    pub expression_pairwise_classifier: Option<f64>,
}

impl ToJson for CuppaPrediction {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("cancerType", &self.cancer_type)
            .f_double("likelihood", self.likelihood)
            .f_opt_double("snvPairwiseClassifier", self.snv_pairwise_classifier)
            .f_opt_double("genomicPositionClassifier", self.genomic_position_classifier)
            .f_opt_double("featureClassifier", self.feature_classifier)
            .f_opt_double("altSjCohortClassifier", self.alt_sj_cohort_classifier)
            .f_opt_double("expressionPairwiseClassifier", self.expression_pairwise_classifier)
            .end_object();
    }
}

#[derive(Clone, Debug)]
pub struct CuppaData {
    pub mode: String,
    pub predictions: Vec<CuppaPrediction>,
    pub best_prediction: CuppaPrediction,
    pub simple_dups_32_to_200b: i32,
    pub max_complex_size: i32,
    pub telomeric_sgls: i32,
    pub line_count: i32,
}

impl ToJson for CuppaData {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("mode", &self.mode)
            .f_list("predictions", &self.predictions)
            .f_value("bestPrediction", &self.best_prediction)
            .f_int("simpleDups32To200B", self.simple_dups_32_to_200b as i64)
            .f_int("maxComplexSize", self.max_complex_size as i64)
            .f_int("telomericSGLs", self.telomeric_sgls as i64)
            .f_int("LINECount", self.line_count as i64)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PeachGenotype {
    pub gene: String,
    pub allele: String,
    pub allele_count: i32,
    pub haplotype: String,
    pub function: String,
    pub linked_drugs: String,
    pub url_prescription_info: String,
    pub panel_version: String,
    pub repo_version: String,
}

impl ToJson for PeachGenotype {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("gene", &self.gene)
            .f_str("allele", &self.allele)
            .f_int("alleleCount", self.allele_count as i64)
            .f_str("haplotype", &self.haplotype)
            .f_str("function", &self.function)
            .f_str("linkedDrugs", &self.linked_drugs)
            .f_str("urlPrescriptionInfo", &self.url_prescription_info)
            .f_str("panelVersion", &self.panel_version)
            .f_str("repoVersion", &self.repo_version)
            .end_object();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SignatureAllocation {
    pub signature: String,
    pub etiology: String,
    pub allocation: f64,
    pub percent: f64,
}

impl ToJson for SignatureAllocation {
    fn write_json(&self, w: &mut JsonWriter) {
        w.begin_object()
            .f_str("signature", &self.signature)
            .f_str("etiology", &self.etiology)
            .f_double("allocation", self.allocation)
            .f_double("percent", self.percent)
            .end_object();
    }
}

/// OrangeJson.write: the JSON text of the record
pub fn to_json(record: &OrangeRecord) -> String {
    let mut w = JsonWriter::new();
    record.write_json(&mut w);
    w.finish()
}
