//! PaveConstants

#[allow(dead_code)]
pub const APP_NAME: &str = "Pave";

pub const GENE_UPSTREAM_DISTANCE: i32 = 1000;
pub const SPLICE_REGION_INTRON_RANGE: i32 = 8;
pub const SPLICE_REGION_EXON_RANGE: i32 = 3;

pub const SPLICE_ACCEPTOR_END_RANGE: i32 = 3;
pub const SPLICE_ACCEPTOR_POSITIONS: [i32; 3] = [1, 2, SPLICE_ACCEPTOR_END_RANGE];
pub const SPLICE_DONOR_POSITIONS: [i32; 4] = [-1, 1, 2, 5];

/// only TERT sets a coding impact upstream of the coding start
pub const PROMOTOR_UPSTREAM_GENE_IDS: [&str; 1] = ["ENSG00000164362"];
pub const PROMOTOR_UPSTREAM_DISTANCE: i32 = 300;

// PON
pub const PON_SAMPLE_COUNT_THRESHOLD: i32 = 10;
pub const PON_REPEAT_COUNT_THRESHOLD: i32 = 4;
pub const PON_MEAN_READ_THRESHOLD: i32 = 6;
pub const PON_VAF_THRESHOLD: f64 = 0.08;

pub const PON_INDEL_ARTEFACT_REPEAT_COUNT: i32 = 7;
pub const PON_INDEL_ARTEFACT_VAF_REDUCTION: f64 = 0.1;
pub const PON_INDEL_ARTEFACT_SBX_FACTOR: i32 = 40;
pub const PON_INDEL_ARTEFACT_SBX_MIN_THRESHOLD: i32 = 8;
pub const PON_INDEL_ARTEFACT_ULTIMA_FACTOR: i32 = 60;
pub const PON_INDEL_ARTEFACT_ULTIMA_MIN_THRESHOLD: i32 = 18;

pub const GNMOAD_FILTER_HOTSPOT_PATHOGENIC_THRESHOLD: f64 = 0.01;
pub const GNMOAD_FILTER_THRESHOLD: f64 = 0.00015;

// filters and VCF tags
pub const PASS_FILTER: &str = "PASS";
pub const PON_FILTER: &str = "PON";
pub const PON_ARTEFACT_FILTER: &str = "PONArtefact";
pub const PON_GNOMAD_FILTER: &str = "PONGnomad";

pub const PON_COUNT: &str = "PON_COUNT";
pub const PON_MAX: &str = "PON_MAX";
pub const PON_AVG_READS: &str = "PON_READS";
pub const PON_MULTI_STATUS: &str = "PON_MULTI";

pub const GNOMAD_FREQ: &str = "GND_FREQ";
pub const GNOMAD_FREQ_DESC: &str = "Gnomad variant frequency";
pub const MAPPABILITY: &str = "MAPPABILITY";
pub const MAPPABILITY_DESC: &str = "GEM mappability in 150 base window";
pub const BLACKLIST_BED_FLAG: &str = "BLACKLIST_BED";
pub const BLACKLIST_VCF_FLAG: &str = "BLACKLIST_VCF";
pub const BLACKLIST_VCF_DESC: &str = "Sites listed in Blacklist VCF";
pub const CLNSIG: &str = "CLNSIG";
pub const CLNSIGCONF: &str = "CLNSIGCONF";
pub const CLNSIG_DESC: &str = "Clinical significance for this single variant";
pub const CLNSIGCONF_DESC: &str = "Conflicting clinical significance for this single variant";
pub const HOTSPOT_FLAG: &str = "HOTSPOT";
pub const HOTSPOT_DESCRIPTION: &str = "Site is at a known hotspot location";
pub const REPORTED_FLAG: &str = "REPORTED";
pub const REPORTED_DESC: &str = "Variant is reported in the driver catalog";

pub const VAR_TRANS_IMPACT_ANNOTATION: &str = "PAVE_TI";
pub const VAR_IMPACT: &str = "IMPACT";
pub const PROTEIN_CONTEXT: &str = "PR_CXT";
pub const PROTEIN_CONTEXT_DESC: &str = "Protein context information";

/// PonCache.PON_FILTERS_V38 (the 1000-sample PON) and PON_FILTERS_V37
pub const PON_FILTERS_V38: &str = "HOTSPOT:6:5;PANEL:3:3;UNKNOWN:3:0";
pub const PON_FILTERS_V37: &str = "HOTSPOT:10:5;PANEL:6:5;UNKNOWN:6:0";
