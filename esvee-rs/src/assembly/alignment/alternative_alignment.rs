//! AlternativeAlignment: one of the alternative locations BWA reports in the XA tag.

use crate::assembly::common_utils::human_chromosome_contains;
use crate::assembly::orientation::Orientation;
use crate::assembly::sam::get_alignment_end_from_cigar;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlternativeAlignment {
    pub chromosome: String,
    pub position: i32,
    pub orient: Orientation,
    pub cigar: String,
    pub num_mutations: i32,
}

/// String.split(delim, limit) with a positive limit
fn split_limit(s: &str, delim: char, limit: usize) -> Vec<&str> {
    s.splitn(limit, delim).collect()
}

impl AlternativeAlignment {
    pub fn new(chromosome: &str, position: i32, orientation: Orientation, cigar: &str, num_mutations: i32) -> Self {
        AlternativeAlignment { chromosome: chromosome.to_string(), position, orient: orientation, cigar: cigar.to_string(), num_mutations }
    }
    pub fn position_end(&self) -> i32 {
        get_alignment_end_from_cigar(self.position, &self.cigar)
    }

    fn from_bwa_string(mapping_str: &str) -> Option<AlternativeAlignment> {
        let items = split_limit(mapping_str, ',', 4);
        if items.len() != 4 {
            return None;
        }
        let chromosome = items[0];
        if !human_chromosome_contains(chromosome) {
            return None;
        }
        let orient_pos = items[1];
        let orientation = Orientation::from_char(orient_pos.chars().next().expect("XA orientation"));
        let position: i32 = orient_pos[1..].parse().expect("XA position");
        let num_mutations: i32 = items[3].parse().expect("XA edit distance");
        Some(AlternativeAlignment::new(chromosome, position, orientation, items[2], num_mutations))
    }

    pub fn from_location_tag(location_tag: &str) -> Vec<AlternativeAlignment> {
        if location_tag.is_empty() {
            return Vec::new();
        }
        location_tag.split(';').filter_map(Self::from_bwa_string).collect()
    }

    pub fn vcf_string(&self) -> String {
        format!("{}:{}|{}|{}|{}", self.chromosome, self.position, self.orient.as_char(), self.cigar, self.num_mutations)
    }

    pub fn to_vcf_tag(alts: &[AlternativeAlignment]) -> String {
        alts.iter().map(|x| x.vcf_string()).collect::<Vec<_>>().join(",")
    }

    pub fn from_vcf_tag(vcf_tag: &str) -> Vec<AlternativeAlignment> {
        // StructuralVariantFactory.trimStringListValue: strips enclosing square brackets
        let trimmed = vcf_tag.strip_prefix('[').unwrap_or(vcf_tag);
        let trimmed = trimmed.strip_suffix(']').unwrap_or(trimmed);
        if trimmed.is_empty() {
            return Vec::new();
        }
        trimmed.split(',').filter_map(Self::from_vcf_string).collect()
    }

    fn from_vcf_string(alignment_str: &str) -> Option<AlternativeAlignment> {
        let values: Vec<&str> = alignment_str.split('|').collect();
        if values.len() != 4 {
            return None;
        }
        let mut location: Vec<&str> = values[0].split(':').collect();
        while location.last().is_some_and(|x| x.is_empty()) {
            location.pop();
        }
        if location.len() != 2 {
            return None;
        }
        if values[1].is_empty() {
            return None;
        }
        let chromosome = location[0];
        if !human_chromosome_contains(chromosome) {
            return None;
        }
        let position: i32 = location[1].parse().expect("alt alignment position");
        let orientation = Orientation::from_char(values[1].chars().next().unwrap());
        let num_mutations: i32 = values[3].parse().expect("alt alignment edit distance");
        Some(AlternativeAlignment::new(chromosome, position, orientation, values[2], num_mutations))
    }

    pub fn matches(&self, other: &AlternativeAlignment) -> bool {
        self == other
    }
}

impl std::fmt::Display for AlternativeAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:{}:{} nm={}", self.chromosome, self.position, self.orient.as_byte(), self.num_mutations)
    }
}
