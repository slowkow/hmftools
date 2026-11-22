//! Junction: a candidate breakend position and orientation from the prep step.

use std::cmp::Ordering;
use std::collections::HashMap;

use anyhow::{Context, Result};

use crate::assembly::common_utils::{compare_junctions, IndelCoords};
use crate::assembly::constants::MIN_VARIANT_LENGTH;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::{A, P};
use crate::assembly::sam::intern;
use crate::common::java_sort::java_sort;
use crate::common::regions::SpecificRegions;

pub type JuncP = P<Junction>;

pub const INDEL_TYPE: &str = "I";
pub const HOTSPOT_TYPE: &str = "H";
pub const DISCORDANT_TYPE: &str = "D";

pub struct Junction {
    pub chromosome: &'static str,
    pub position: i32,
    pub orient: Orientation,
    pub discordant_only: bool,
    pub hotspot: bool,
    pub max_soft_clip_length: i32,
    indel_based: bool,
    pub details: String,
    raw_discordant_position: i32,
    indel_coords: Option<A<IndelCoords>>,
    saga_match_id: Option<String>,
}

impl Junction {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        chromosome: &str, position: i32, orientation: Orientation, discordant_only: bool, indel_based: bool, hotspot: bool,
        max_soft_clip_length: i32, saga_match_id: Option<String>,
    ) -> Junction {
        let details = if discordant_only || indel_based || hotspot {
            let mut parts = Vec::new();
            if discordant_only {
                parts.push("disc-only");
            }
            if indel_based {
                parts.push("indel");
            }
            if hotspot {
                parts.push("hotspot");
            }
            parts.join("/")
        } else {
            String::new()
        };
        Junction {
            chromosome: intern(chromosome),
            position,
            orient: orientation,
            discordant_only,
            hotspot,
            max_soft_clip_length,
            indel_based,
            details,
            raw_discordant_position: 0,
            indel_coords: None,
            saga_match_id,
        }
    }

    /// the test constructor Junction(chromosome, position, orientation)
    pub fn new_basic(chromosome: &str, position: i32, orientation: Orientation) -> Junction {
        let mut j = Junction::new(chromosome, position, orientation, false, false, false, MIN_VARIANT_LENGTH, None);
        j.raw_discordant_position = -1;
        j
    }

    pub fn is_forward(&self) -> bool {
        self.orient.is_forward()
    }
    pub fn is_reverse(&self) -> bool {
        self.orient.is_reverse()
    }
    pub fn set_indel_status(&mut self, is_indel: bool) {
        self.indel_based = is_indel;
    }
    pub fn indel_based(&self) -> bool {
        self.indel_based
    }
    pub fn soft_clip_based(&self) -> bool {
        !self.indel_based() && !self.discordant_only
    }
    pub fn set_indel_coords(&mut self, indel_coords: Option<A<IndelCoords>>) {
        self.indel_coords = indel_coords;
    }
    pub fn indel_coords(&self) -> Option<&A<IndelCoords>> {
        self.indel_coords.as_ref()
    }
    pub fn raw_discordant_position(&self) -> i32 {
        self.raw_discordant_position
    }
    pub fn set_raw_discordant_position(&mut self, position: i32) {
        self.raw_discordant_position = position;
    }

    pub fn coords(&self) -> String {
        format!("{}:{}:{}", self.chromosome, self.position, self.orient.as_byte())
    }
    pub fn coords_typed(&self) -> String {
        self.coords_typed_ext(false)
    }
    pub fn coords_typed_ext(&self, use_raw_discordant: bool) -> String {
        if self.discordant_only {
            format!("{}:{}:{}:D", self.chromosome, if use_raw_discordant { self.raw_discordant_position } else { self.position }, self.orient.as_byte())
        } else if self.indel_based() {
            format!("{}:I", self.coords())
        } else {
            self.coords()
        }
    }

    pub fn is_local_match(&self, other: &Junction) -> bool {
        self.position == other.position && self.orient == other.orient
    }
    pub fn is_saga_matched(&self) -> bool {
        self.saga_match_id.is_some()
    }
    pub fn saga_match(&self) -> Option<&str> {
        self.saga_match_id.as_deref()
    }
    pub fn lower(&self, other: &Junction) -> bool {
        self.position < other.position || (self.position == other.position && self.orient.is_forward() && other.orient.is_reverse())
    }
    pub fn higher(&self, other: &Junction) -> bool {
        !self.lower(other)
    }
    pub fn compare_to(&self, other: &Junction) -> Ordering {
        compare_junctions(self.chromosome, other.chromosome, self.position, other.position, self.orient, other.orient)
    }

    /// Junction.loadJunctions. The map is keyed by chromosome; each list is sorted.
    pub fn load_junctions(
        filename: &str, specific_regions: &SpecificRegions, min_junction_frags: i32, min_hotspot_frags: i32, min_discordant_frags: i32,
    ) -> Result<HashMap<String, Vec<JuncP>>> {
        let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read junctions file({filename})"))?;
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let field = |name: &str| header.iter().position(|h| *h == name);
        let chr_index = field("Chromosome").context("junction file: no Chromosome column")?;
        let pos_index = field("Position").context("junction file: no Position column")?;
        let orient_index = field("Orientation").context("junction file: no Orientation column")?;
        let junc_frags_index = field("JunctionFrags");
        let other_junc_frags_index = field("ExactSupportFrags");
        let other_support_frags_index = field("OtherSupportFrags").or_else(|| field("DiscordantFrags"));
        let indel_index = field("Indel");
        let hotspot_index = field("Hotspot");
        let extra_info_index = field("ExtraInfo");
        let remote_frags_index = field("RemoteFrags");
        let saga_match_index = field("SagaMatch");
        let max_sc_length_index = field("MaxValidSoftClip");

        let mut chr_junctions_map: HashMap<String, Vec<JuncP>> = HashMap::new();
        let mut current_chromosome = String::new();
        let mut junction_count = 0;
        let mut raw_junction_count = 0;
        let mut discordant_count = 0;
        let mut indel_count = 0;
        let mut hotspot_count = 0;
        let int = |s: &str| -> Result<i32> { s.parse::<i32>().with_context(|| format!("junction file: invalid number '{s}'")) };
        let parse_bool = |s: &str| s.eq_ignore_ascii_case("true");

        for line in lines {
            let values: Vec<&str> = line.split('\t').collect();
            let chromosome = values[chr_index];
            let position = int(values[pos_index])?;
            let orientation = Orientation::from_byte_str(values[orient_index]);
            if specific_regions.has_filters() {
                if !specific_regions.include_chromosome(chromosome) {
                    continue;
                }
                if !specific_regions.regions.is_empty() && !specific_regions.regions.iter().any(|r| r.chromosome == chromosome && r.contains_position(position)) {
                    continue;
                }
            }
            raw_junction_count += 1;
            let junction_frags = match junc_frags_index {
                Some(i) => int(values[i])?,
                None => 0,
            };
            // as in Java, the presence of the other-support column decides whether the exact-support column is read
            let other_junction_frags = match other_support_frags_index {
                Some(_) => int(values[other_junc_frags_index.context("junction file: no ExactSupportFrags column")?])?,
                None => 0,
            };
            let other_support_frags = match other_support_frags_index {
                Some(i) => int(values[i])?,
                None => 0,
            };
            let max_soft_clip_length = match max_sc_length_index {
                Some(i) => int(values[i])?,
                None => MIN_VARIANT_LENGTH,
            };
            let discordant_only = junction_frags == 0 && other_support_frags > 0;
            let indel = indel_index.is_some_and(|i| parse_bool(values[i]));
            let hotspot = hotspot_index.is_some_and(|i| parse_bool(values[i]));
            let saga_match_variant = saga_match_index.map(|i| values[i].split(' ').next().unwrap_or("").to_string());

            if hotspot {
                if junction_frags < min_hotspot_frags {
                    continue;
                }
                hotspot_count += 1;
            } else if discordant_only {
                let mut max_remote_frags = other_support_frags;
                if let Some(i) = remote_frags_index {
                    max_remote_frags = int(values[i])?;
                } else if let Some(i) = extra_info_index {
                    max_remote_frags = int(values[i])?;
                }
                if max_remote_frags < min_discordant_frags {
                    continue;
                }
                discordant_count += 1;
            } else {
                if junction_frags + other_junction_frags < min_junction_frags {
                    continue;
                }
                if indel {
                    indel_count += 1;
                }
            }
            if current_chromosome != chromosome {
                current_chromosome = chromosome.to_string();
                chr_junctions_map.insert(current_chromosome.clone(), Vec::new());
            }
            chr_junctions_map.get_mut(&current_chromosome).unwrap().push(P::new(Junction::new(
                chromosome, position, orientation, discordant_only, indel, hotspot, max_soft_clip_length, saga_match_variant,
            )));
            junction_count += 1;
        }
        let split_count = junction_count - discordant_count - indel_count;
        crate::info!(
            "loaded junctions({} raw={}) types(split={} indel={} discordant={} hotspot={}) from file: {}",
            junction_count, raw_junction_count, split_count, indel_count, discordant_count, hotspot_count, filename
        );
        for list in chr_junctions_map.values_mut() {
            java_sort(list, |a, b| a.compare_to(b));
        }
        Ok(chr_junctions_map)
    }

    pub fn from_config_str(config_str: &str) -> Option<Junction> {
        let values: Vec<&str> = config_str.split(':').collect();
        if values.len() < 3 {
            return None;
        }
        let chromosome = values[0];
        let position: i32 = values[1].parse().unwrap_or_else(|_| panic!("invalid junction position: {}", values[1]));
        let orientation = Orientation::from_byte_str(values[2]);
        let mut discordant_only = false;
        let mut indel_based = false;
        let mut hotspot = false;
        if values.len() >= 4 {
            if values[3] == DISCORDANT_TYPE {
                discordant_only = true;
            } else if values[3] == INDEL_TYPE {
                indel_based = true;
            }
        }
        if values.len() >= 5 && values[4] == HOTSPOT_TYPE {
            hotspot = true;
        }
        Some(Junction::new(chromosome, position, orientation, discordant_only, indel_based, hotspot, MIN_VARIANT_LENGTH, None))
    }

    /// Junction.mergeJunctions (used by tests)
    pub fn merge_junctions(existing_map: &mut HashMap<String, Vec<JuncP>>, new_map: &HashMap<String, Vec<JuncP>>) {
        if existing_map.is_empty() {
            for (k, v) in new_map {
                existing_map.insert(k.clone(), v.clone());
            }
            return;
        }
        for (new_chromosome, new_junctions) in new_map {
            if new_junctions.is_empty() {
                continue;
            }
            let existing_junctions = match existing_map.get_mut(new_chromosome) {
                Some(e) if !e.is_empty() => e,
                _ => {
                    existing_map.insert(new_chromosome.clone(), new_junctions.clone());
                    continue;
                }
            };
            let mut existing_index = 0usize;
            let mut new_index = 0usize;
            let mut existing_junction: Option<JuncP> = Some(existing_junctions[existing_index]);
            let mut new_junction = new_junctions[new_index];
            loop {
                let mut move_new = false;
                let mut move_existing = false;
                match existing_junction {
                    None => {
                        existing_junctions.push(new_junction);
                        move_new = true;
                    }
                    Some(ej) => {
                        if ej.is_local_match(&new_junction) {
                            move_new = true;
                            move_existing = true;
                        } else if new_junction.higher(&ej) {
                            move_existing = true;
                        } else if new_junction.lower(&ej) {
                            existing_junctions.insert(existing_index, new_junction);
                            existing_index += 1;
                            move_new = true;
                        } else {
                            crate::error!("junction merge failed: existing({}:{}) vs new({}:{})", existing_index, *ej, new_index, *new_junction);
                            return;
                        }
                    }
                }
                if move_new {
                    new_index += 1;
                    if new_index >= new_junctions.len() {
                        break;
                    }
                    new_junction = new_junctions[new_index];
                }
                if move_existing {
                    existing_index += 1;
                    existing_junction = existing_junctions.get(existing_index).copied();
                }
            }
        }
    }

    pub fn validate_junction_map(map: &HashMap<String, Vec<JuncP>>) -> bool {
        for junctions in map.values() {
            for i in 0..junctions.len().saturating_sub(1) {
                let junction = junctions[i];
                let next_junction = junctions[i + 1];
                if junction.is_local_match(&next_junction) {
                    crate::error!("junction({}:{}) duplicate", i, *junction);
                    return false;
                }
                if junction.position > next_junction.position {
                    crate::error!("junction({}:{}) after next({})", i, *junction, *next_junction);
                    return false;
                }
                for (j, other) in junctions.iter().enumerate().skip(i + 1) {
                    if junction.is_local_match(other) {
                        crate::error!("junction({}:{}) duplicate with later index({})", i, *junction, j);
                        return false;
                    }
                }
            }
        }
        true
    }
}

impl std::fmt::Display for Junction {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.discordant_only || self.indel_based || self.hotspot {
            write!(f, "{}:{}:{} {}", self.chromosome, self.position, self.orient.as_byte(), self.details)
        } else {
            write!(f, "{}:{}:{}", self.chromosome, self.position, self.orient.as_byte())
        }
    }
}
