//! Panel-of-normals matching (PonCache, PonSvRegion, PonSglRegion).

use super::filters::*;
use super::types::*;
use super::variant::*;
use anyhow::{bail, Context, Result};
use hashbrown::HashMap;
use std::io::BufRead;

#[derive(Clone, Debug, PartialEq)]
pub struct PonSvRegion {
    pub start: i32,
    pub end: i32,
    pub orient_start: Orient,
    /// index into PonCache::end_chromosomes
    pub end_chromosome: u16,
    pub end_start: i32,
    pub end_end: i32,
    pub orient_end: Orient,
    pub pon_count: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PonSglRegion {
    pub start: i32,
    pub end: i32,
    pub orient: Orient,
    pub pon_count: i32,
}

pub struct PonCache {
    /// SV entries by start chromosome, SGL entries by chromosome
    sv_regions: HashMap<String, Vec<PonSvRegion>>,
    sgl_regions: HashMap<String, Vec<PonSglRegion>>,
    /// end chromosome names without a `chr` prefix, as the end regions compare chromosomes by number
    end_chromosomes: Vec<String>,
    position_sv_margin: i32,
    position_sgl_margin: i32,
    current_sv_index: usize,
    current_sgl_index: usize,
    has_valid_data: bool,
    specific_sgl_fusion_regions: Vec<ChrBaseRegion>,
    /// log lines of the loading, issued by the caller in order
    pub messages: Vec<(bool, String)>,
}

fn open_lines(path: &str) -> Result<Box<dyn BufRead>> {
    let file = super::vcf::open_text(path).with_context(|| path.to_string())?;
    Ok(Box::new(std::io::BufReader::with_capacity(1 << 18, file)))
}

fn parse_int(s: &str) -> Result<i32> {
    super::vcf::java_parse_int(s).with_context(|| format!("NumberFormatException: For input string: \"{s}\""))
}

fn orient_from_char(s: &str) -> Result<Orient> {
    match s.chars().next() {
        Some('+') => Ok(Orient::Forward),
        Some('-') => Ok(Orient::Reverse),
        other => bail!("Invalid orientation: {other:?}"),
    }
}

fn orient_from_byte_str(s: &str) -> Result<Orient> {
    match s {
        "1" | "+1" => Ok(Orient::Forward),
        "-1" => Ok(Orient::Reverse),
        other => bail!("Invalid orientation: {other}"),
    }
}

impl PonCache {
    /// PonCache(svMargin, sglMargin, ponSvFile, ponSglFile, allowUnordered = false)
    pub fn new(sv_margin: i32, sgl_margin: i32, sv_file: Option<&str>, sgl_file: Option<&str>) -> Result<PonCache> {
        let mut cache = PonCache {
            sv_regions: HashMap::new(),
            sgl_regions: HashMap::new(),
            end_chromosomes: Vec::new(),
            position_sv_margin: sv_margin,
            position_sgl_margin: sgl_margin,
            current_sv_index: 0,
            current_sgl_index: 0,
            has_valid_data: true,
            specific_sgl_fusion_regions: Vec::new(),
            messages: Vec::new(),
        };
        if let Some(file) = sv_file {
            cache.load_sv_file(file)?;
        }
        if let Some(file) = sgl_file {
            cache.load_sgl_file(file)?;
        }
        Ok(cache)
    }

    /// regions where a SGL with an alternative alignment is not matched to the PON
    pub fn set_specific_sgl_regions(&mut self, version: RefGenomeVersion) {
        self.specific_sgl_fusion_regions = if version.is37() {
            vec![
                ChrBaseRegion::new("4", 190930000, 191030000),
                ChrBaseRegion::new("10", 135420000, 135520000),
                ChrBaseRegion::new("X", 52729628, 52731680),
                ChrBaseRegion::new("X", 52784877, 52786929),
            ]
        } else {
            vec![
                ChrBaseRegion::new("chr4", 190060000, 190190000),
                ChrBaseRegion::new("chr10", 133660000, 133770000),
                ChrBaseRegion::new("chrX", 52700578, 52702630),
                ChrBaseRegion::new("chrX", 52755800, 52757852),
            ]
        };
    }

    pub fn has_valid_data(&self) -> bool {
        self.has_valid_data
    }

    fn end_chromosome_index(&mut self, chromosome: &str) -> Result<u16> {
        let name = human_chromosome_name(chromosome);
        if let Some(i) = self.end_chromosomes.iter().position(|c| *c == name) {
            return Ok(i as u16);
        }
        if self.end_chromosomes.len() >= u16::MAX as usize {
            bail!("too many distinct end chromosomes in SV PON");
        }
        self.end_chromosomes.push(name);
        Ok((self.end_chromosomes.len() - 1) as u16)
    }

    #[cfg(test)]
    pub fn add_sv_region(&mut self, chr_start: &str, start: (i32, i32), orient_start: Orient, chr_end: &str, end: (i32, i32), orient_end: Orient, pon_count: i32) {
        let end_chromosome = self.end_chromosome_index(chr_end).unwrap();
        self.sv_regions.entry(chr_start.to_string()).or_default().push(PonSvRegion {
            start: start.0,
            end: start.1,
            orient_start,
            end_chromosome,
            end_start: end.0,
            end_end: end.1,
            orient_end,
            pon_count,
        });
    }

    #[cfg(test)]
    pub fn add_sgl_region(&mut self, chromosome: &str, region: (i32, i32), orient: Orient, pon_count: i32) {
        self.sgl_regions.entry(chromosome.to_string()).or_default().push(PonSglRegion { start: region.0, end: region.1, orient, pon_count });
    }

    fn load_sv_file(&mut self, filename: &str) -> Result<()> {
        let reader = match open_lines(filename) {
            Ok(r) => r,
            Err(e) => {
                self.messages.push((true, format!("failed to load SV PON file({filename}): {e:#}")));
                self.has_valid_data = false;
                return Ok(());
            }
        };
        let is_tsv = filename.ends_with(".tsv.gz");
        let mut item_count = 0;
        let mut current_chr = String::new();
        let mut regions: Vec<PonSvRegion> = Vec::new();
        let mut last_start: Option<i32> = None;
        let mut end_index_cache: (String, u16) = (String::new(), 0);

        for (n, line) in reader.lines().enumerate() {
            let line = match line {
                Ok(l) => l,
                Err(e) => {
                    self.messages.push((true, format!("failed to load SV PON file({filename}): {e}")));
                    self.has_valid_data = false;
                    return Ok(());
                }
            };
            if is_tsv && n == 0 {
                continue;
            }
            let items: Vec<&str> = line.split('\t').collect();
            let field = |i: usize| items.get(i).copied().with_context(|| format!("SV PON line with {} fields: {line}", items.len()));
            let (chr_start, start, end, chr_end, end_start, end_end, orient_start, orient_end, pon_count);
            if is_tsv {
                chr_start = field(0)?;
                start = parse_int(field(1)?)?;
                end = parse_int(field(2)?)?;
                chr_end = field(3)?;
                end_start = parse_int(field(4)?)?;
                end_end = parse_int(field(5)?)?;
                orient_start = orient_from_byte_str(field(6)?)?;
                orient_end = orient_from_byte_str(field(7)?)?;
                pon_count = parse_int(field(8)?)?;
            } else {
                chr_start = field(0)?;
                start = parse_int(field(1)?)? + 1;
                end = parse_int(field(2)?)?;
                chr_end = field(3)?;
                end_start = parse_int(field(4)?)? + 1;
                end_end = parse_int(field(5)?)?;
                orient_start = orient_from_char(field(8)?)?;
                orient_end = orient_from_char(field(9)?)?;
                pon_count = parse_int(field(7)?)?;
            }

            if chr_start != current_chr {
                if !current_chr.is_empty() {
                    self.sv_regions.entry(std::mem::take(&mut current_chr)).or_default().append(&mut regions);
                }
                current_chr = chr_start.to_string();
                last_start = None;
            }
            if end_index_cache.0 != chr_end {
                end_index_cache = (chr_end.to_string(), self.end_chromosome_index(chr_end)?);
            }
            regions.push(PonSvRegion { start, end, orient_start, end_chromosome: end_index_cache.1, end_start, end_end, orient_end, pon_count });
            item_count += 1;

            if let Some(last) = last_start {
                if last > start {
                    let previous = &regions[regions.len() - 2];
                    self.messages.push((
                        false,
                        format!("SV PON not ordered: last({current_chr}:{}-{}) vs this({current_chr}:{start}-{end})", previous.start, previous.end),
                    ));
                    self.has_valid_data = false;
                }
            }
            last_start = Some(start);
        }
        if !current_chr.is_empty() {
            self.sv_regions.entry(current_chr).or_default().append(&mut regions);
        }
        for list in self.sv_regions.values_mut() {
            list.shrink_to_fit();
        }
        self.messages.push((false, format!("loaded {item_count} SV PON records from file({filename})")));
        Ok(())
    }

    fn load_sgl_file(&mut self, filename: &str) -> Result<()> {
        let reader = match open_lines(filename) {
            Ok(r) => r,
            Err(e) => {
                self.messages.push((true, format!("failed to load SGL PON file({filename}): {e:#}")));
                self.has_valid_data = false;
                return Ok(());
            }
        };
        let is_tsv = filename.ends_with(".tsv.gz");
        let mut item_count = 0;
        let mut current_chr = String::new();
        let mut regions: Vec<PonSglRegion> = Vec::new();
        let mut last_start: Option<i32> = None;

        for (n, line) in reader.lines().enumerate() {
            let line = match line {
                Ok(l) => l,
                Err(e) => {
                    self.messages.push((true, format!("failed to load SGL PON file({filename}): {e}")));
                    self.has_valid_data = false;
                    return Ok(());
                }
            };
            if is_tsv && n == 0 {
                continue;
            }
            let items: Vec<&str> = line.split('\t').collect();
            let field = |i: usize| items.get(i).copied().with_context(|| format!("SGL PON line with {} fields: {line}", items.len()));
            let (chromosome, start, end, orient, pon_count);
            if is_tsv {
                chromosome = field(0)?;
                start = parse_int(field(1)?)?;
                end = parse_int(field(2)?)?;
                orient = orient_from_byte_str(field(3)?)?;
                pon_count = parse_int(field(4)?)?;
            } else {
                chromosome = field(0)?;
                start = parse_int(field(1)?)? + 1;
                end = parse_int(field(2)?)?;
                orient = orient_from_char(field(5)?)?;
                pon_count = parse_int(field(4)?)?;
            }

            if chromosome != current_chr {
                if !current_chr.is_empty() {
                    self.sgl_regions.entry(std::mem::take(&mut current_chr)).or_default().append(&mut regions);
                }
                current_chr = chromosome.to_string();
                last_start = None;
            }
            regions.push(PonSglRegion { start, end, orient, pon_count });
            item_count += 1;

            if let Some(last) = last_start {
                if last > start {
                    let previous = &regions[regions.len() - 2];
                    self.messages.push((
                        false,
                        format!("SGL PON not ordered: last({current_chr}:{}-{}) vs this({current_chr}:{start}-{end})", previous.start, previous.end),
                    ));
                    self.has_valid_data = false;
                }
            }
            last_start = Some(start);
        }
        if !current_chr.is_empty() {
            self.sgl_regions.entry(current_chr).or_default().append(&mut regions);
        }
        for list in self.sgl_regions.values_mut() {
            list.shrink_to_fit();
        }
        self.messages.push((false, format!("loaded {item_count} SGL PON records from file({filename})")));
        Ok(())
    }

    /// PonCache.annotateVariants
    pub fn annotate_variants(&mut self, variants: &mut [Variant], map: &BreakendMap) {
        for (chromosome, breakends) in map.iter() {
            self.current_sv_index = 0;
            self.current_sgl_index = 0;
            let sgl_regions = self.sgl_regions.remove(chromosome.as_str());
            let sv_regions = self.sv_regions.remove(chromosome.as_str());

            for b in breakends {
                let var = &variants[b.var as usize];
                // only looked up on the first breakend
                if !var.is_sgl() && b.se as usize == SE_END {
                    continue;
                }
                if var.pon_count > 0 {
                    continue;
                }
                let pon_count =
                    if var.is_sgl() { self.find_sgl_pon_match(sgl_regions.as_deref(), var) } else { self.find_sv_pon_match(sv_regions.as_deref(), var) };
                if pon_count > 0 {
                    let var = &mut variants[b.var as usize];
                    var.pon_count = pon_count;
                    var.add_filter(FilterType::Pon);
                }
            }
            if let Some(regions) = sgl_regions {
                self.sgl_regions.insert(chromosome.clone(), regions);
            }
            if let Some(regions) = sv_regions {
                self.sv_regions.insert(chromosome.clone(), regions);
            }
        }
    }

    /// the margins either side of a breakend's position: (start, end)
    fn breakend_margin(var: &Variant, se: usize) -> (i32, i32) {
        let breakend = var.breakend(se);
        let inexact_homology = (if breakend.is_start { breakend.inexact_homology.start } else { breakend.inexact_homology.end }).abs();
        let mut margin_start = breakend.confidence_interval.start.min(-inexact_homology);
        let mut margin_end = breakend.confidence_interval.end.max(inexact_homology);
        if !var.is_sgl() {
            let max_ins_seq_length = (breakend.insert_sequence.len() as i32).min(PON_MAX_INS_SEQ_LENGTH);
            if breakend.orient.is_forward() {
                margin_end += max_ins_seq_length;
            } else {
                margin_start -= max_ins_seq_length;
            }
        }
        (margin_start, margin_end)
    }

    fn sv_position_margin(&self, var: &Variant) -> i32 {
        let mut position_margin = self.position_sv_margin;
        if var.is_short_local() && var.sv_length() <= PON_SHORT_INDEL_LENGTH {
            position_margin = PON_SHORT_INDEL_PON_DISTANCE;
            let max_local_repeat = var.start().context.attr_int(MAX_LOCAL_REPEAT, 0);
            if max_local_repeat >= PON_SHORT_INDEL_MAX_REPEATS {
                position_margin = PON_SHORT_INDEL_MAX_REPEAT_PON_DISTANCE;
            }
        }
        position_margin
    }

    fn find_sv_pon_match(&mut self, regions: Option<&[PonSvRegion]>, var: &Variant) -> i32 {
        let Some(regions) = regions else { return 0 };
        let (mut start_lo, mut start_hi) = Self::breakend_margin(var, SE_START);
        let (mut end_lo, mut end_hi) = Self::breakend_margin(var, SE_END);
        let position_margin = self.sv_position_margin(var);
        start_lo -= position_margin;
        start_hi += position_margin;
        end_lo -= position_margin;
        end_hi += position_margin;

        let sv_start = (var.pos_start() + start_lo, var.pos_start() + start_hi);

        while self.current_sv_index < regions.len() {
            let region = &regions[self.current_sv_index];
            if positions_overlap(region.start, region.end, sv_start.0, sv_start.1) {
                let sv_end = (var.pos_end() + end_lo, var.pos_end() + end_hi);
                return self.find_pon_match(regions, var, sv_start, sv_end, self.current_sv_index);
            }
            // exit if the PON is now past this point and retreat one position
            if region.start > sv_start.1 {
                break;
            }
            self.current_sv_index += 1;
        }
        if self.current_sv_index > 0 {
            self.current_sv_index -= 1;
        }
        0
    }

    fn find_pon_match(&self, regions: &[PonSvRegion], var: &Variant, sv_start: (i32, i32), sv_end: (i32, i32), start_index: usize) -> i32 {
        let end_chromosome = human_chromosome_name(var.chromosome_end());
        let end_chromosome_index = self.end_chromosomes.iter().position(|c| *c == end_chromosome);
        let (orient_start, orient_end) = (var.orient_start(), var.orient_end());
        let mut max_pon_count = 0;
        for search_up in [true, false] {
            let mut current = start_index as i64;
            while current >= 0 && (current as usize) < regions.len() {
                let region = &regions[current as usize];
                if search_up && region.start > sv_start.1 {
                    break;
                }
                if !search_up && region.end < sv_start.0 {
                    break;
                }
                if positions_overlap(region.start, region.end, sv_start.0, sv_start.1)
                    && Some(region.end_chromosome as usize) == end_chromosome_index
                    && positions_overlap(region.end_start, region.end_end, sv_end.0, sv_end.1)
                    && region.orient_start == orient_start
                    && region.orient_end == orient_end
                    && region.pon_count > max_pon_count
                {
                    max_pon_count = region.pon_count;
                }
                current += if search_up { 1 } else { -1 };
            }
        }
        max_pon_count
    }

    fn matches_specific_sgl_fusion_region(&self, var: &Variant) -> bool {
        let alignments_str = var.start().context.attr_string(INSALN, "");
        if alignments_str.is_empty() {
            return false;
        }
        AlternativeAlignment::from_vcf_tag(&alignments_str)
            .iter()
            .any(|alt| self.specific_sgl_fusion_regions.iter().any(|r| r.contains_position(&alt.chromosome, alt.position)))
    }

    fn find_sgl_pon_match(&mut self, regions: Option<&[PonSglRegion]>, var: &Variant) -> i32 {
        if self.matches_specific_sgl_fusion_region(var) {
            return 0;
        }
        let Some(regions) = regions else { return 0 };
        let (mut lo, mut hi) = Self::breakend_margin(var, SE_START);
        lo -= self.position_sgl_margin;
        hi += self.position_sgl_margin;
        let sv_start = (var.pos_start() + lo, var.pos_start() + hi);

        while self.current_sgl_index < regions.len() {
            let region = &regions[self.current_sgl_index];
            if positions_overlap(region.start, region.end, sv_start.0, sv_start.1) {
                return Self::find_sgl_match_around(regions, var, sv_start, self.current_sgl_index);
            }
            // exit if the PON is now past this point then retreat one position
            if region.start > sv_start.1 {
                if self.current_sgl_index > 0 {
                    self.current_sgl_index -= 1;
                }
                break;
            }
            self.current_sgl_index += 1;
        }
        0
    }

    fn find_sgl_match_around(regions: &[PonSglRegion], var: &Variant, sv_start: (i32, i32), start_index: usize) -> i32 {
        let orient = var.orient_start();
        let mut max_pon_count = 0;
        for search_up in [true, false] {
            let mut current = start_index as i64;
            while current >= 0 && (current as usize) < regions.len() {
                let region = &regions[current as usize];
                if search_up && region.start > sv_start.1 {
                    break;
                }
                if !search_up && region.end < sv_start.0 {
                    break;
                }
                if positions_overlap(region.start, region.end, sv_start.0, sv_start.1) && region.orient == orient && region.pon_count > max_pon_count {
                    max_pon_count = region.pon_count;
                }
                current += if search_up { 1 } else { -1 };
            }
        }
        max_pon_count
    }
}
