//! hmf-common `fusion`: the known fusion reference data (KnownFusionCache, KnownFusionData,
//! KnownFusionType) and FusionCommon's constants.

use anyhow::{bail, Result};

use crate::java_text::{fields_index_map, parse_int, read_lines};
use crate::regions::ChrBaseRegion;

pub const KNOWN_FUSIONS_FILE: &str = "known_fusion_file";
pub const KNOWN_FUSIONS_FILE_DESC: &str = "Known fusion reference data file";

pub const DEFAULT_PRE_GENE_PROMOTOR_DISTANCE: i32 = 100000;
pub const FS_UP: usize = 0;
pub const FS_DOWN: usize = 1;
pub const FS_PAIR: usize = 2;
pub const POS_STRAND: i8 = 1;
pub const NEG_STRAND: i8 = -1;

pub fn fs_index(is_upstream: bool) -> usize {
    if is_upstream { FS_UP } else { FS_DOWN }
}
pub fn switch_stream(i: usize) -> usize {
    if i == FS_UP { FS_DOWN } else { FS_UP }
}
pub fn stream_str(i: usize) -> &'static str {
    if i == FS_UP { "upstream" } else { "downstream" }
}

/// KnownFusionType
#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum KnownFusionType {
    NONE,
    KNOWN_PAIR,
    PROMISCUOUS_5,
    PROMISCUOUS_3,
    ENHANCER_KNOWN_PAIR,
    ENHANCER_PROMISCUOUS,
    EXON_DEL_DUP,
    PROMISCUOUS_ENHANCER_TARGET,
}

pub const PROMISCUOUS_BOTH: &str = "PROMISCUOUS_BOTH";

impl KnownFusionType {
    pub const VALUES: [KnownFusionType; 8] = [
        Self::NONE,
        Self::KNOWN_PAIR,
        Self::PROMISCUOUS_5,
        Self::PROMISCUOUS_3,
        Self::ENHANCER_KNOWN_PAIR,
        Self::ENHANCER_PROMISCUOUS,
        Self::EXON_DEL_DUP,
        Self::PROMISCUOUS_ENHANCER_TARGET,
    ];
    pub fn name(self) -> String {
        format!("{self:?}")
    }
    pub fn value_of(s: &str) -> Option<Self> {
        Self::VALUES.iter().copied().find(|t| t.name() == s)
    }
    /// KnownFusionType.parse: the pre-3.0 IG names too
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "IG_KNOWN_PAIR" => Some(Self::ENHANCER_KNOWN_PAIR),
            "IG_PROMISCUOUS" => Some(Self::ENHANCER_PROMISCUOUS),
            _ => Self::value_of(s),
        }
    }
}

impl std::fmt::Display for KnownFusionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// KnownFusionData
#[derive(Clone, Debug)]
pub struct KnownFusionData {
    pub fusion_type: KnownFusionType,
    pub five_gene: String,
    pub three_gene: String,
    pub cancer_types: String,
    pub pub_med_id: String,
    valid_data: bool,
    high_impact_promiscuous: bool,
    downstream_distance: [i32; 2],
    specific_exons_trans_name: String,
    five_exon_range: [i32; 2],
    three_exon_range: [i32; 2],
    gene_region: Option<ChrBaseRegion>,
    /// None: either orientation
    gene_orientation: Option<i8>,
    three_gene_alt_regions: Vec<ChrBaseRegion>,
}

impl KnownFusionData {
    pub fn new(t: KnownFusionType, five_gene: &str, three_gene: &str, cancer_types: &str, pub_med_id: &str) -> Self {
        KnownFusionData {
            fusion_type: t,
            five_gene: five_gene.to_string(),
            three_gene: three_gene.to_string(),
            cancer_types: cancer_types.to_string(),
            pub_med_id: pub_med_id.to_string(),
            valid_data: true,
            high_impact_promiscuous: false,
            downstream_distance: [0, 0],
            specific_exons_trans_name: String::new(),
            five_exon_range: [0, 0],
            three_exon_range: [0, 0],
            gene_region: None,
            gene_orientation: None,
            three_gene_alt_regions: Vec::new(),
        }
    }

    /// KnownFusionData.fromCsv: an error where Java's constructor arguments throw (the caller skips the
    /// line); errors in the type-specific fields are logged and leave the rest set
    pub fn from_csv(data: &str, f: &std::collections::HashMap<String, usize>) -> Result<KnownFusionData> {
        let items: Vec<&str> = data.split(',').collect();
        let get = |n: &str| -> Result<&str> {
            let Some(i) = f.get(n) else { bail!("missing column {n}") };
            items.get(*i).copied().ok_or_else(|| anyhow::anyhow!("ArrayIndexOutOfBounds"))
        };
        let t = KnownFusionType::value_of(get("Type")?).ok_or_else(|| anyhow::anyhow!("invalid known fusion type"))?;
        let mut kf = KnownFusionData::new(t, get("FiveGene")?, get("ThreeGene")?, get("CancerTypes")?, get("PubMedId")?);
        let known_exon_trans = get("KnownExonTranscript")?;
        let up = get("KnownExonUpRange")?;
        let down = get("KnownExonDownRange")?;
        let r = (|| -> Result<()> {
            if !known_exon_trans.is_empty() {
                kf.set_known_exon_data(known_exon_trans, up, down)?;
            }
            if get("HighImpactPromiscuous")?.eq_ignore_ascii_case("TRUE") {
                kf.high_impact_promiscuous = true;
            }
            let overrides = get("Overrides")?;
            if !overrides.is_empty() {
                kf.apply_overrides(overrides)?;
            }
            Ok(())
        })();
        if let Err(e) = r {
            crate::error!("failed to parse specific data for known fusion({}): error({})", kf, e);
        }
        Ok(kf)
    }

    pub fn set_known_exon_data(&mut self, trans: &str, up: &str, down: &str) -> Result<()> {
        self.specific_exons_trans_name = trans.to_string();
        let range = |s: &str| -> Result<[i32; 2]> {
            let e = crate::java_text::split(s, ';');
            let a = e.first().and_then(|x| parse_int(x)).ok_or_else(|| anyhow::anyhow!("NumberFormatException: {s}"))?;
            let b = e.get(1).and_then(|x| parse_int(x)).ok_or_else(|| anyhow::anyhow!("ArrayIndexOutOfBounds: {s}"))?;
            Ok([a, b])
        };
        if !up.is_empty() {
            self.five_exon_range = range(up)?;
        }
        if !down.is_empty() {
            self.three_exon_range = range(down)?;
        }
        Ok(())
    }

    pub fn apply_overrides(&mut self, overrides: &str) -> Result<()> {
        for item in crate::java_text::split(overrides, ' ') {
            let parts = crate::java_text::split(item, '=');
            let name = parts.first().copied().unwrap_or("");
            let Some(data) = parts.get(1).copied() else { bail!("ArrayIndexOutOfBounds: {item}") };
            let int = |s: &str| parse_int(s).ok_or_else(|| anyhow::anyhow!("NumberFormatException: For input string: \"{s}\""));
            match name {
                "UP_GENE_DOWNSTREAM_DISTANCE" => self.downstream_distance[FS_UP] = int(data)?,
                "DOWN_GENE_DOWNSTREAM_DISTANCE" => self.downstream_distance[FS_DOWN] = int(data)?,
                "ALTS" => {
                    let alt = crate::java_text::split(data, ';');
                    let mut index = 0;
                    while index < alt.len() {
                        if alt[index] != "ALT" || alt.len() - index < 4 {
                            self.valid_data = false;
                            return Ok(());
                        }
                        index += 1;
                        self.three_gene_alt_regions.push(ChrBaseRegion::new(alt[index], int(alt[index + 1])?, int(alt[index + 2])?));
                        index += 3;
                    }
                }
                "ENHANCER_RANGE" | "THREE_PRIME_RANGE" => {
                    let r: Vec<&str> = data.split(';').collect();
                    let o = r.first().copied().unwrap_or("");
                    self.gene_orientation = if o == "0" {
                        None
                    } else {
                        Some(match o {
                            "1" => 1,
                            "-1" => -1,
                            _ => bail!("invalid orientation: {o}"),
                        })
                    };
                    if r.len() < 4 {
                        bail!("ArrayIndexOutOfBounds: {data}");
                    }
                    self.gene_region = Some(ChrBaseRegion::new(r[1], int(r[2])?, int(r[3])?));
                }
                _ => {}
            }
        }
        Ok(())
    }

    pub fn valid_data(&self) -> bool {
        self.valid_data
    }
    pub fn is_high_impact_promiscuous(&self) -> bool {
        self.high_impact_promiscuous
    }
    pub fn downstream_distance(&self, stream: usize) -> i32 {
        self.downstream_distance[stream]
    }
    pub fn specific_exons_trans_name(&self) -> &str {
        &self.specific_exons_trans_name
    }
    pub fn five_gene_exon_range(&self) -> [i32; 2] {
        self.five_exon_range
    }
    pub fn three_gene_exon_range(&self) -> [i32; 2] {
        self.three_exon_range
    }
    pub fn gene_region(&self) -> Option<&ChrBaseRegion> {
        self.gene_region.as_ref()
    }
    pub fn gene_strand(&self) -> Option<i8> {
        self.gene_orientation
    }
    pub fn within_gene_region(&self, chromosome: &str, position: i32) -> bool {
        self.gene_region.as_ref().is_some_and(|r| r.chromosome == chromosome && r.contains_position(position))
    }
    pub fn three_gene_alt_regions(&self) -> &[ChrBaseRegion] {
        &self.three_gene_alt_regions
    }
}

impl std::fmt::Display for KnownFusionData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let range = self.gene_region.as_ref().map(|r| r.to_string()).unwrap_or_else(|| "none".into());
        write!(f, "{}: genes({} - {}) range({}) ct({})", self.fusion_type, self.five_gene, self.three_gene, range, self.cancer_types)
    }
}

/// KnownFusionCache
#[derive(Clone)]
pub struct KnownFusionCache {
    data: Vec<KnownFusionData>,
    has_valid_data: bool,
}

impl Default for KnownFusionCache {
    fn default() -> Self {
        Self::new()
    }
}

impl KnownFusionCache {
    pub fn new() -> Self {
        KnownFusionCache { data: Vec::new(), has_valid_data: true }
    }
    pub fn has_valid_data(&self) -> bool {
        self.has_valid_data
    }
    pub fn data(&self) -> &[KnownFusionData] {
        &self.data
    }
    pub fn data_by_type(&self, t: KnownFusionType) -> impl Iterator<Item = &KnownFusionData> {
        self.data.iter().filter(move |x| x.fusion_type == t)
    }
    fn known_pairs(&self) -> impl Iterator<Item = &KnownFusionData> {
        self.data_by_type(KnownFusionType::KNOWN_PAIR)
    }
    pub fn add_data(&mut self, d: KnownFusionData) {
        self.data.push(d);
    }

    pub fn has_known_fusion(&self, five: &str, three: &str) -> bool {
        self.known_pairs().any(|x| x.five_gene == five && x.three_gene == three)
    }
    pub fn has_known_unmappable3_fusion(&self, five: &str, three: &str) -> bool {
        self.known_pairs().filter(|x| !x.three_gene_alt_regions.is_empty()).any(|x| x.five_gene == five && x.three_gene == three)
    }
    pub fn is_single_breakend_candidate(&self, gene: &str, is_upstream: bool) -> bool {
        if self.known_pairs().any(|x| (is_upstream && x.five_gene == gene) || (!is_upstream && x.three_gene == gene)) {
            return true;
        }
        self.data_by_type(KnownFusionType::ENHANCER_KNOWN_PAIR)
            .filter(|x| !x.three_gene_alt_regions.is_empty())
            .any(|x| !is_upstream && x.three_gene == gene)
    }
    pub fn has_promiscuous_five_gene(&self, gene: &str) -> bool {
        self.data_by_type(KnownFusionType::PROMISCUOUS_5).any(|x| x.five_gene == gene)
    }
    pub fn has_promiscuous_three_gene(&self, gene: &str) -> bool {
        self.data_by_type(KnownFusionType::PROMISCUOUS_3).any(|x| x.three_gene == gene)
    }
    pub fn has_any_enhancer_fusion(&self, gene: &str) -> bool {
        self.data_by_type(KnownFusionType::ENHANCER_KNOWN_PAIR).any(|x| x.five_gene == gene || x.three_gene == gene)
    }
    pub fn has_exon_del_dup(&self, gene: &str) -> bool {
        self.data_by_type(KnownFusionType::EXON_DEL_DUP).any(|x| x.five_gene == gene && x.three_gene == gene)
    }
    pub fn has_known_pair_gene(&self, gene: &str) -> bool {
        self.known_pairs().any(|x| x.five_gene == gene || x.three_gene == gene)
    }
    pub fn is_exon_del_dup_trans(&self, trans_name: &str) -> bool {
        self.data_by_type(KnownFusionType::EXON_DEL_DUP).any(|x| x.specific_exons_trans_name == trans_name)
    }
    pub fn requires_promiscuous_exon_range(&self, t: KnownFusionType, trans_name: &str) -> bool {
        for k in self.data_by_type(t) {
            if k.specific_exons_trans_name != trans_name {
                continue;
            }
            let r = if t == KnownFusionType::PROMISCUOUS_5 { k.five_exon_range } else { k.three_exon_range };
            if r[0] >= 1 {
                return true;
            }
        }
        false
    }
    pub fn within_promiscuous_exon_range(&self, t: KnownFusionType, trans_name: &str, breakend_exon: i32, fused_exon: i32) -> bool {
        for k in self.data_by_type(t) {
            if k.specific_exons_trans_name != trans_name {
                continue;
            }
            let r = if t == KnownFusionType::PROMISCUOUS_5 { k.five_exon_range } else { k.three_exon_range };
            if breakend_exon >= r[0] && breakend_exon <= r[1] && fused_exon >= r[0] && fused_exon <= r[1] {
                return true;
            }
        }
        false
    }
    pub fn within_known_exon_ranges(
        &self, t: KnownFusionType, trans_name: &str, breakend_exon_up: i32, fused_exon_up: i32, breakend_exon_down: i32, fused_exon_down: i32,
    ) -> bool {
        for k in self.data_by_type(t) {
            if k.specific_exons_trans_name != trans_name {
                continue;
            }
            let f = k.five_exon_range;
            let th = k.three_exon_range;
            if breakend_exon_up >= f[0]
                && breakend_exon_up <= f[1]
                && fused_exon_up >= f[0]
                && fused_exon_up <= f[1]
                && breakend_exon_down >= th[0]
                && breakend_exon_down <= th[1]
                && fused_exon_down >= th[0]
                && fused_exon_down <= th[1]
            {
                return true;
            }
        }
        false
    }
    pub fn is_high_impact_promiscuous(&self, t: KnownFusionType, five: &str, three: &str) -> bool {
        for k in self.data.iter().filter(|x| x.high_impact_promiscuous) {
            if k.fusion_type != t {
                continue;
            }
            if k.fusion_type == KnownFusionType::PROMISCUOUS_5 && k.five_gene == five {
                return true;
            }
            if k.fusion_type == KnownFusionType::PROMISCUOUS_3 && k.three_gene == three {
                return true;
            }
        }
        false
    }
    pub fn within_enhancer_region(&self, chromosome: &str, position: i32) -> bool {
        self.data.iter().filter(|x| x.gene_region.is_some()).any(|x| x.within_gene_region(chromosome, position))
    }

    /// loadFromFile: false (and the cache invalid) if the file is missing or unreadable
    pub fn load_from_file(&mut self, filename: &str) -> bool {
        if !self.load_file(filename) {
            self.has_valid_data = false;
            return false;
        }
        let counts: Vec<String> = KnownFusionType::VALUES
            .iter()
            .filter(|t| **t != KnownFusionType::NONE)
            .map(|t| format!("{}({})", t, self.data_by_type(*t).count()))
            .collect();
        crate::info!("loaded known fusion data: {}", counts.join(", "));
        true
    }

    fn load_file(&mut self, filename: &str) -> bool {
        if !std::path::Path::new(filename).exists() {
            crate::error!("file({}) not found", filename);
            return false;
        }
        let lines = match read_lines(filename) {
            Ok(l) => l,
            Err(e) => {
                crate::error!("file({}) invalid known fusion data: {}", filename, e);
                return false;
            }
        };
        if lines.is_empty() {
            return false;
        }
        let f = fields_index_map(&lines[0], ',');
        for line in lines.iter().skip(1) {
            match KnownFusionData::from_csv(line, &f) {
                Ok(d) => self.add_data(d),
                Err(_) => crate::warn!("file({}) invalid known fusion data will be skipped: {}", filename, line),
            }
        }
        true
    }
}
