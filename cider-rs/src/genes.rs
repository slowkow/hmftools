//! The IG/TCR genes: loci, V/J gene types, gene locations (genes/*.kt), the gene file of hmf-common
//! (IgTcrGeneFile, igtcr_gene.<37|38>.tsv), CiderGeneDataLoader and CiderGeneDatastore.

use std::collections::HashMap;

use anyhow::{anyhow, bail, Result};
use hmf_common::regions::RefGenomeVersion;

use crate::compat;

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum Vj {
    V,
    J,
}

/// IgTcrLocus, in declaration order
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum IgTcrLocus {
    Igh,
    Igk,
    Igl,
    TraTrd,
    Trb,
    Trg,
}

impl IgTcrLocus {
    pub const ALL: [IgTcrLocus; 6] = [IgTcrLocus::Igh, IgTcrLocus::Igk, IgTcrLocus::Igl, IgTcrLocus::TraTrd, IgTcrLocus::Trb, IgTcrLocus::Trg];

    pub fn pretty_print(self) -> &'static str {
        match self {
            IgTcrLocus::Igh => "IGH",
            IgTcrLocus::Igk => "IGK",
            IgTcrLocus::Igl => "IGL",
            IgTcrLocus::TraTrd => "TRA/TRD",
            IgTcrLocus::Trb => "TRB",
            IgTcrLocus::Trg => "TRG",
        }
    }

    /// fromGeneName: the first 3 characters; TRA and TRD are one locus
    pub fn from_gene_name(gene_name: &str) -> Result<IgTcrLocus> {
        let name = gene_name.get(0..3).ok_or_else(|| anyhow!("invalid gene name: {gene_name}"))?;
        Ok(match name {
            "TRA" | "TRD" => IgTcrLocus::TraTrd,
            "IGH" => IgTcrLocus::Igh,
            "IGK" => IgTcrLocus::Igk,
            "IGL" => IgTcrLocus::Igl,
            "TRB" => IgTcrLocus::Trb,
            "TRG" => IgTcrLocus::Trg,
            _ => bail!("No enum constant com.hartwig.hmftools.cider.genes.IgTcrLocus.{name}"),
        })
    }
}

/// VJGeneType, in declaration order (the order of EnumMap and of sorting)
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum VjGeneType {
    Ighv,
    Ighj,
    Igkv,
    Igkj,
    Iglv,
    Iglj,
    Trav,
    Traj,
    Trbv,
    Trbj,
    Trdv,
    Trdj,
    Trgv,
    Trgj,
}

impl VjGeneType {
    pub const ALL: [VjGeneType; 14] = [
        VjGeneType::Ighv,
        VjGeneType::Ighj,
        VjGeneType::Igkv,
        VjGeneType::Igkj,
        VjGeneType::Iglv,
        VjGeneType::Iglj,
        VjGeneType::Trav,
        VjGeneType::Traj,
        VjGeneType::Trbv,
        VjGeneType::Trbj,
        VjGeneType::Trdv,
        VjGeneType::Trdj,
        VjGeneType::Trgv,
        VjGeneType::Trgj,
    ];

    pub fn ordinal(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            VjGeneType::Ighv => "IGHV",
            VjGeneType::Ighj => "IGHJ",
            VjGeneType::Igkv => "IGKV",
            VjGeneType::Igkj => "IGKJ",
            VjGeneType::Iglv => "IGLV",
            VjGeneType::Iglj => "IGLJ",
            VjGeneType::Trav => "TRAV",
            VjGeneType::Traj => "TRAJ",
            VjGeneType::Trbv => "TRBV",
            VjGeneType::Trbj => "TRBJ",
            VjGeneType::Trdv => "TRDV",
            VjGeneType::Trdj => "TRDJ",
            VjGeneType::Trgv => "TRGV",
            VjGeneType::Trgj => "TRGJ",
        }
    }

    pub fn locus(self) -> IgTcrLocus {
        use VjGeneType::*;
        match self {
            Ighv | Ighj => IgTcrLocus::Igh,
            Igkv | Igkj => IgTcrLocus::Igk,
            Iglv | Iglj => IgTcrLocus::Igl,
            Trav | Traj | Trdv | Trdj => IgTcrLocus::TraTrd,
            Trbv | Trbj => IgTcrLocus::Trb,
            Trgv | Trgj => IgTcrLocus::Trg,
        }
    }

    pub fn vj(self) -> Vj {
        if self.ordinal() % 2 == 0 { Vj::V } else { Vj::J }
    }

    /// pairedVjGeneTypes
    pub fn paired(self) -> &'static [VjGeneType] {
        use VjGeneType::*;
        match self {
            Ighv => &[Ighj],
            Ighj => &[Ighv],
            Igkv => &[Igkj],
            Igkj => &[Igkv],
            Iglv => &[Iglj],
            Iglj => &[Iglv],
            Trav => &[Traj, Trdj],
            Traj => &[Trav, Trdv],
            Trbv => &[Trbj],
            Trbj => &[Trbv],
            Trdv => &[Trdj, Traj],
            Trdj => &[Trdv, Trav],
            Trgv => &[Trgj],
            Trgj => &[Trgv],
        }
    }

    /// fromGeneName: IGKINTR is IGKV, IGKDEL is IGKJ, else the first 4 characters
    pub fn from_gene_name(gene_name: &str) -> Result<VjGeneType> {
        match gene_name {
            "IGKINTR" => return Ok(VjGeneType::Igkv),
            "IGKDEL" => return Ok(VjGeneType::Igkj),
            _ => {}
        }
        let prefix = gene_name.get(0..4).ok_or_else(|| anyhow!("invalid gene name: {gene_name}"))?;
        VjGeneType::ALL
            .iter()
            .copied()
            .find(|t| t.name() == prefix)
            .ok_or_else(|| anyhow!("No enum constant com.hartwig.hmftools.cider.genes.VJGeneType.{prefix}"))
    }
}

/// hmf-common IgTcrRegion
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum IgTcrRegion {
    VRegion,
    DRegion,
    JRegion,
    Constant,
}

impl IgTcrRegion {
    pub fn is_vj(self) -> bool {
        matches!(self, IgTcrRegion::VRegion | IgTcrRegion::JRegion)
    }
    pub fn is_vdj(self) -> bool {
        self.is_vj() || self == IgTcrRegion::DRegion
    }
    fn from_name(s: &str) -> Result<IgTcrRegion> {
        Ok(match s {
            "V_REGION" => IgTcrRegion::VRegion,
            "D_REGION" => IgTcrRegion::DRegion,
            "J_REGION" => IgTcrRegion::JRegion,
            "CONSTANT" => IgTcrRegion::Constant,
            _ => bail!("No enum constant com.hartwig.hmftools.common.cider.IgTcrRegion.{s}"),
        })
    }
}

/// hmf-common IgTcrFunctionality
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum IgTcrFunctionality {
    Functional,
    Orf,
    Pseudogene,
}

impl IgTcrFunctionality {
    pub fn to_code(self) -> &'static str {
        match self {
            IgTcrFunctionality::Functional => "F",
            IgTcrFunctionality::Orf => "ORF",
            IgTcrFunctionality::Pseudogene => "P",
        }
    }
    fn from_code(s: &str) -> Result<IgTcrFunctionality> {
        Ok(match s {
            "F" => IgTcrFunctionality::Functional,
            "ORF" => IgTcrFunctionality::Orf,
            "P" => IgTcrFunctionality::Pseudogene,
            _ => bail!("invalid IgTcrFunctionality code: {s}"),
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn as_char(self) -> char {
        if self == Strand::Forward { '+' } else { '-' }
    }
}

/// GenomicLocation: 1-based inclusive positions; possibly in a non-primary assembly
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GenomicLocation {
    pub chromosome: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub strand: Strand,
    pub alt_assembly_name: Option<String>,
}

impl GenomicLocation {
#[allow(dead_code)]
    pub fn new(chromosome: &str, pos_start: i32, pos_end: i32, strand: Strand) -> GenomicLocation {
        assert!(pos_start <= pos_end);
        GenomicLocation { chromosome: chromosome.to_string(), pos_start, pos_end, strand, alt_assembly_name: None }
    }
    pub fn in_primary_assembly(&self) -> bool {
        self.alt_assembly_name.is_none()
    }
    pub fn base_length(&self) -> i32 {
        self.pos_end - self.pos_start + 1
    }
}

impl std::fmt::Display for GenomicLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(a) = &self.alt_assembly_name {
            write!(f, "{a} ")?;
        }
        write!(f, "{}:{}-{}({})", self.chromosome, self.pos_start, self.pos_end, self.strand.as_char())
    }
}

/// IgTcrGene (cider), from a row of the gene file
#[derive(Clone, Debug)]
pub struct IgTcrGene {
    pub gene_name: String,
    pub allele: String,
    pub region: IgTcrRegion,
    pub functionality: IgTcrFunctionality,
    pub gene_location: Option<GenomicLocation>,
    pub anchor_sequence: Option<String>,
    pub anchor_location: Option<GenomicLocation>,
}

impl IgTcrGene {
    pub fn gene_allele(&self) -> String {
        format!("{}*{}", self.gene_name, self.allele)
    }
    pub fn is_functional(&self) -> bool {
        self.functionality == IgTcrFunctionality::Functional
    }
}

pub const GENE_FILE_38: &str = include_str!("../resources/igtcr_gene.38.tsv");
pub const GENE_FILE_37: &str = include_str!("../resources/igtcr_gene.37.tsv");

/// IgTcrGeneFile.read, then IgTcrGene.fromCommonIgTcrGene.
///
/// CI1: fromCommonIgTcrGene builds the gene and anchor locations from the region and strand only, without
/// the assembly name, so the locations of genes on ALT, patch and unplaced contigs (whose coordinates
/// are positions within those contigs, written against the primary chromosome's name) count as primary
/// assembly. Corrected: the assembly name is kept, so those locations are left out of the read
/// selection, as CiderGeneDataLoader intends.
pub fn read_gene_file(version: RefGenomeVersion) -> Result<Vec<IgTcrGene>> {
    let text = if version.is_37() { GENE_FILE_37 } else { GENE_FILE_38 };
    parse_gene_file(text, version, !compat::java("CI1"))
}

pub fn parse_gene_file(text: &str, version: RefGenomeVersion, keep_assembly_name: bool) -> Result<Vec<IgTcrGene>> {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().ok_or_else(|| anyhow!("empty gene file"))?.split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).ok_or_else(|| anyhow!("missing column {name}"));
    let (c_gene, c_allele, c_region, c_func, c_primary, c_asm, c_chr, c_start, c_end, c_strand, c_astart, c_aend, c_aseq) = (
        col("gene")?,
        col("allele")?,
        col("region")?,
        col("functionality")?,
        col("primaryAssembly")?,
        col("assemblyName")?,
        col("chromosome")?,
        col("posStart")?,
        col("posEnd")?,
        col("strand")?,
        col("anchorStart")?,
        col("anchorEnd")?,
        col("anchorSequence")?,
    );
    let mut genes = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let region = IgTcrRegion::from_name(get(c_region))?;
        let functionality = IgTcrFunctionality::from_code(get(c_func))?;
        // getBooleanOrNull: empty is null (primary)
        let primary = match get(c_primary) {
            "" => None,
            s => Some(s.eq_ignore_ascii_case("true")),
        };
        let assembly_name = if primary.is_none() || primary == Some(true) { None } else { Some(get(c_asm).to_string()) };
        let anchor_sequence = if get(c_aseq).is_empty() { None } else { Some(get(c_aseq).to_string()) };
        let mut chromosome = get(c_chr).to_string();
        let mut gene_location = None;
        let mut anchor_location = None;
        if !chromosome.is_empty() {
            chromosome = version.versioned_chromosome(&chromosome);
            let pos_start: i32 = get(c_start).parse()?;
            let pos_end: i32 = get(c_end).parse()?;
            let strand = match get(c_strand).chars().next() {
                Some('+') => Strand::Forward,
                Some('-') => Strand::Reverse,
                None => bail!("chromosome exists but strand is missing"),
                Some(c) => bail!("invalid strand: {c}"),
            };
            if pos_start <= 0 || pos_end <= 0 {
                bail!("chromosome exists but pos start or pos end is invalid");
            }
            let alt = if keep_assembly_name { assembly_name.clone() } else { None };
            gene_location = Some(GenomicLocation { chromosome: chromosome.clone(), pos_start, pos_end, strand, alt_assembly_name: alt.clone() });
            let (a_start, a_end) = (get(c_astart), get(c_aend));
            if !a_start.is_empty() && !a_end.is_empty() {
                anchor_location =
                    Some(GenomicLocation { chromosome: chromosome.clone(), pos_start: a_start.parse()?, pos_end: a_end.parse()?, strand, alt_assembly_name: alt });
            }
        }
        genes.push(IgTcrGene {
            gene_name: get(c_gene).to_string(),
            allele: get(c_allele).to_string(),
            region,
            functionality,
            gene_location,
            anchor_sequence,
            anchor_location,
        });
    }
    Ok(genes)
}

/// VJAnchorTemplate
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct VjAnchorTemplate {
    pub gene_type: VjGeneType,
    pub gene_name: String,
    pub allele: String,
    pub gene_location: Option<GenomicLocation>,
    pub anchor_sequence: String,
    pub anchor_location: Option<GenomicLocation>,
    pub anchor_amino_acid_sequence: String,
}

impl VjAnchorTemplate {
    pub fn new(
        gene_type: VjGeneType,
        gene_name: &str,
        allele: &str,
        gene_location: Option<GenomicLocation>,
        anchor_sequence: &str,
        anchor_location: Option<GenomicLocation>,
    ) -> VjAnchorTemplate {
        assert!(!anchor_sequence.is_empty());
        assert!(anchor_location.as_ref().is_none_or(|l| l.in_primary_assembly()));
        VjAnchorTemplate {
            gene_type,
            gene_name: gene_name.to_string(),
            allele: allele.to_string(),
            gene_location,
            anchor_sequence: anchor_sequence.to_string(),
            anchor_location,
            anchor_amino_acid_sequence: hmf_common::codon::amino_acid_from_bases(anchor_sequence.as_bytes()),
        }
    }
}

/// VJAnchorGenomeLocation
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct VjAnchorGenomeLocation {
    pub gene_type: VjGeneType,
    pub location: GenomicLocation,
}

impl VjAnchorGenomeLocation {
    pub fn vj(&self) -> Vj {
        self.gene_type.vj()
    }
    pub fn start(&self) -> i32 {
        self.location.pos_start
    }
    pub fn end(&self) -> i32 {
        self.location.pos_end
    }
    pub fn strand(&self) -> Strand {
        self.location.strand
    }
    pub fn base_length(&self) -> i32 {
        self.location.base_length()
    }
    /// the side of the anchor towards the CDR3: 1 for the end of the C codon of V on the forward strand
    /// and the start of J on the reverse strand, else -1
    pub fn anchor_boundary_side(&self) -> i32 {
        let right = (self.vj() == Vj::V && self.strand() == Strand::Forward) || (self.vj() == Vj::J && self.strand() == Strand::Reverse);
        if right { 1 } else { -1 }
    }
    pub fn anchor_boundary_reference_position(&self) -> i32 {
        if self.anchor_boundary_side() == 1 { self.end() } else { self.start() }
    }
}

impl std::fmt::Display for VjAnchorGenomeLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "VJAnchorGenomeLocation(vjGeneType={}, genomeLocation={})", self.gene_type.name(), self.location)
    }
}

/// IgTcrConstantDiversityRegion
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct IgTcrConstantDiversityRegion {
    pub locus: IgTcrLocus,
    pub location: GenomicLocation,
    pub gene_name: String,
}

impl IgTcrConstantDiversityRegion {
    pub fn corresponding_vj(&self) -> &'static [VjGeneType] {
        use VjGeneType::*;
        match self.locus {
            IgTcrLocus::Igh => &[Ighv, Ighj],
            IgTcrLocus::Igk => &[Igkv, Igkj],
            IgTcrLocus::Igl => &[Iglv, Iglj],
            IgTcrLocus::Trb => &[Trbv, Trbj],
            IgTcrLocus::Trg => &[Trgv, Trgj],
            IgTcrLocus::TraTrd => &[Trav, Trdv, Traj, Trdj],
        }
    }
}

/// CiderGeneDataLoader.loadConstantDiversityRegions
pub fn load_constant_diversity_regions(genes: &[IgTcrGene]) -> Result<Vec<IgTcrConstantDiversityRegion>> {
    let mut regions: Vec<IgTcrConstantDiversityRegion> = Vec::new();
    for g in genes {
        let Some(loc) = &g.gene_location else { continue };
        if !loc.in_primary_assembly() {
            continue;
        }
        if g.region != IgTcrRegion::Constant && g.region != IgTcrRegion::DRegion {
            continue;
        }
        let r = IgTcrConstantDiversityRegion { locus: IgTcrLocus::from_gene_name(&g.gene_name)?, location: loc.clone(), gene_name: g.gene_name.clone() };
        if regions.contains(&r) {
            continue;
        }
        regions.push(r);
    }
    Ok(regions)
}

fn primary_assembly_location(l: &Option<GenomicLocation>) -> Option<GenomicLocation> {
    l.as_ref().filter(|l| l.in_primary_assembly()).cloned()
}

/// CiderGeneDataLoader.loadAnchorTemplates
pub fn load_anchor_templates(genes: &[IgTcrGene]) -> Result<Vec<VjAnchorTemplate>> {
    let mut templates = Vec::new();
    for g in genes {
        if g.region != IgTcrRegion::VRegion && g.region != IgTcrRegion::JRegion {
            continue;
        }
        let Some(anchor_sequence) = &g.anchor_sequence else { continue };
        let gene_type = VjGeneType::from_gene_name(&g.gene_name)?;
        templates.push(VjAnchorTemplate::new(
            gene_type,
            &g.gene_name,
            &g.allele,
            primary_assembly_location(&g.gene_location),
            anchor_sequence,
            primary_assembly_location(&g.anchor_location),
        ));
    }
    Ok(templates)
}

/// CiderGeneDatastore: templates by anchor sequence and by anchor location, the anchor locations (one per
/// location, in gene file order) and the constant and D regions.
#[allow(dead_code)]
pub struct GeneDatastore {
    pub templates: Vec<VjAnchorTemplate>,
    /// template indices by anchor sequence
    by_anchor_sequence: HashMap<String, Vec<usize>>,
    /// per gene type: template indices by anchor sequence
    by_type_anchor_sequence: Vec<HashMap<String, Vec<usize>>>,
    /// per gene type: the anchor sequences, sorted
    type_anchor_sequences: Vec<Vec<String>>,
    by_location: HashMap<GenomicLocation, Vec<usize>>,
    pub anchor_locations: Vec<VjAnchorGenomeLocation>,
    pub constant_diversity_regions: Vec<IgTcrConstantDiversityRegion>,
}

impl GeneDatastore {
    pub fn new(templates: Vec<VjAnchorTemplate>, constant_diversity_regions: Vec<IgTcrConstantDiversityRegion>) -> Result<GeneDatastore> {
        let mut by_anchor_sequence: HashMap<String, Vec<usize>> = HashMap::new();
        let mut by_type_anchor_sequence: Vec<HashMap<String, Vec<usize>>> = vec![HashMap::new(); VjGeneType::ALL.len()];
        let mut by_location: HashMap<GenomicLocation, Vec<usize>> = HashMap::new();
        let mut location_type: HashMap<GenomicLocation, VjGeneType> = HashMap::new();
        let mut anchor_locations = Vec::new();
        for (i, t) in templates.iter().enumerate() {
            if let Some(loc) = &t.anchor_location {
                by_location.entry(loc.clone()).or_default().push(i);
                match location_type.get(loc) {
                    None => {
                        location_type.insert(loc.clone(), t.gene_type);
                        anchor_locations.push(VjAnchorGenomeLocation { gene_type: t.gene_type, location: loc.clone() });
                    }
                    Some(&existing) if existing != t.gene_type => {
                        crate::error!("gene location: {} is used by multiple gene type: {} and {}", loc, existing.name(), t.gene_type.name());
                        bail!("gene location: {} is used by multiple gene type: {} and {}", loc, existing.name(), t.gene_type.name());
                    }
                    _ => {}
                }
            }
            if !t.anchor_sequence.is_empty() {
                by_anchor_sequence.entry(t.anchor_sequence.clone()).or_default().push(i);
                by_type_anchor_sequence[t.gene_type.ordinal()].entry(t.anchor_sequence.clone()).or_default().push(i);
            }
        }
        let type_anchor_sequences = by_type_anchor_sequence
            .iter()
            .map(|m| {
                let mut v: Vec<String> = m.keys().cloned().collect();
                v.sort();
                v
            })
            .collect();
        crate::info!("found {} gene locations", by_location.len());
        Ok(GeneDatastore { templates, by_anchor_sequence, by_type_anchor_sequence, type_anchor_sequences, by_location, anchor_locations, constant_diversity_regions })
    }

    /// getAnchorSequenceSet(geneType): sorted
    pub fn anchor_sequences(&self, gene_type: VjGeneType) -> &[String] {
        &self.type_anchor_sequences[gene_type.ordinal()]
    }

#[allow(dead_code)]
    pub fn by_anchor_sequence(&self, seq: &str) -> &[usize] {
        self.by_anchor_sequence.get(seq).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn by_type_anchor_sequence(&self, gene_type: VjGeneType, seq: &str) -> &[usize] {
        self.by_type_anchor_sequence[gene_type.ordinal()].get(seq).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn by_location(&self, loc: &GenomicLocation) -> &[usize] {
        self.by_location.get(loc).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn template(&self, i: usize) -> &VjAnchorTemplate {
        &self.templates[i]
    }

    /// CiderGeneDatastore from the gene file of `version`
    pub fn load(version: RefGenomeVersion) -> Result<GeneDatastore> {
        let genes = read_gene_file(version)?;
        let templates = load_anchor_templates(&genes)?;
        let regions = load_constant_diversity_regions(&genes)?;
        GeneDatastore::new(templates, regions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gene_types() {
        assert_eq!(VjGeneType::from_gene_name("IGHV1-2").unwrap(), VjGeneType::Ighv);
        assert_eq!(VjGeneType::from_gene_name("IGKDEL").unwrap(), VjGeneType::Igkj);
        assert_eq!(VjGeneType::from_gene_name("IGKINTR").unwrap(), VjGeneType::Igkv);
        assert_eq!(VjGeneType::Trdj.vj(), Vj::J);
        assert_eq!(VjGeneType::Trgv.vj(), Vj::V);
        assert_eq!(IgTcrLocus::from_gene_name("TRDV1").unwrap(), IgTcrLocus::TraTrd);
        assert!(VjGeneType::from_gene_name("ABCD1").is_err());
    }

    #[test]
    fn gene_file_alt_locations() {
        // CI1: with the assembly name kept, ALT locations are not primary
        let fixed = parse_gene_file(GENE_FILE_38, RefGenomeVersion::V38, true).unwrap();
        let java = parse_gene_file(GENE_FILE_38, RefGenomeVersion::V38, false).unwrap();
        assert_eq!(fixed.len(), java.len());
        let alt_fixed = fixed.iter().filter(|g| g.gene_location.as_ref().is_some_and(|l| !l.in_primary_assembly())).count();
        let alt_java = java.iter().filter(|g| g.gene_location.as_ref().is_some_and(|l| !l.in_primary_assembly())).count();
        assert!(alt_fixed > 0);
        assert_eq!(alt_java, 0);
    }

    #[test]
    fn datastore_loads() {
        for v in [RefGenomeVersion::V37, RefGenomeVersion::V38] {
            let genes = parse_gene_file(if v.is_37() { GENE_FILE_37 } else { GENE_FILE_38 }, v, true).unwrap();
            let ds = GeneDatastore::new(load_anchor_templates(&genes).unwrap(), load_constant_diversity_regions(&genes).unwrap()).unwrap();
            assert!(!ds.anchor_locations.is_empty());
            for t in VjGeneType::ALL {
                let seqs = ds.anchor_sequences(t);
                assert!(seqs.windows(2).all(|w| w[0] < w[1]));
            }
        }
    }
}
