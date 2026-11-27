//! SvUtilities and the small Java helpers the LINX code uses everywhere (Strings.appendStr,
//! CigarUtils.calcCigarAlignedLength, Math.round / min / max, Doubles.round, the reference genome
//! version LINX keeps in a static).

use std::sync::atomic::{AtomicBool, Ordering};

use hmf_common::java_hash::JavaHashMap;

use crate::types::*;

static REF_GENOME_V37: AtomicBool = AtomicBool::new(true);

/// LinxConfig.REF_GENOME_VERSION (V37 until the config sets it)
pub fn set_ref_genome_v37(v37: bool) {
    REF_GENOME_V37.store(v37, Ordering::Relaxed);
}
pub fn ref_genome_is_v37() -> bool {
    REF_GENOME_V37.load(Ordering::Relaxed)
}

/// Strings.appendStr
pub fn append_str(dest: &mut String, source: &str, delim: char) {
    if !dest.is_empty() {
        dest.push(delim);
    }
    dest.push_str(source);
}
pub fn append_string(dest: &str, source: &str, delim: char) -> String {
    if dest.is_empty() { source.to_string() } else { format!("{dest}{delim}{source}") }
}

/// Math.round(double) as a long
pub fn jround(d: f64) -> i64 {
    hmf_common::java_math::java_round(d)
}
/// Math.min(double, double)
pub fn jmin(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a <= b { a } else { b }
}
/// Math.max(double, double)
pub fn jmax(a: f64, b: f64) -> f64 {
    hmf_common::java_math::java_max(a, b)
}

/// Doubles.round(value, decimalPlaces): Math.round(value * 10^n) / 10^n
pub fn doubles_round(value: f64, places: i32) -> f64 {
    let log_value = 10f64.powi(places);
    jround(value * log_value) as f64 / log_value
}

/// Doubles.equal: |a - b| < 1e-10
pub fn doubles_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-10
}

/// CigarUtils.calcCigarAlignedLength: the sum of M, D and N lengths
pub fn calc_cigar_aligned_length(cigar: &str) -> i32 {
    let mut base_length = 0;
    let mut current = 0;
    for c in cigar.chars() {
        if c == 'D' || c == 'M' || c == 'N' {
            base_length += current;
            current = 0;
            continue;
        }
        if let Some(d) = c.to_digit(10) {
            current = current * 10 + d as i32;
        } else {
            current = 0;
        }
    }
    base_length
}

// RefGenomeCoordinates: centromere positions (resources refgenome/centromeres.37.tsv and .38.tsv)
const CENTROMERES_37: [i32; 24] = [
    123035434, 93826171, 92004854, 51160117, 47905641, 60330166, 59554331, 45338887, 48867679, 40754935, 53144205, 36356694, 17500000,
    17500000, 18500000, 36835801, 23763006, 16960898, 26181782, 27869569, 12788129, 14500000, 60132012, 11604553,
];
const CENTROMERES_38: [i32; 24] = [
    123605523, 93139351, 92214016, 50726026, 48272854, 59191911, 59498944, 44955505, 44377363, 40640102, 52751711, 35977330, 17025624,
    17086762, 18362627, 37295920, 24849830, 18161053, 25844927, 28237290, 11890184, 14004553, 60509061, 10430492,
];

/// HumanChromosome index (0 for chr1 .. 23 for Y) of a human chromosome name with or without "chr"
pub fn human_chromosome_index(chromosome: &str) -> Option<usize> {
    if !hmf_common::sv::is_human_chromosome(chromosome) {
        return None;
    }
    let c = hmf_common::regions::strip_chr_prefix(chromosome);
    match c {
        "X" => Some(22),
        "Y" => Some(23),
        _ => c.parse::<usize>().ok().map(|v| v - 1),
    }
}

/// RefGenomeCoordinates.length: 0 for other contigs
pub fn get_chromosome_length(chromosome: &str) -> i32 {
    let Some(_) = human_chromosome_index(chromosome) else { return 0 };
    let c = hmf_common::regions::strip_chr_prefix(chromosome);
    hmf_common::ref_coords::chromosome_length(c, ref_genome_is_v37()).unwrap_or(0)
}

/// RefGenomeCoordinates.centromere: 0 for other contigs
pub fn get_centromere(chromosome: &str) -> i32 {
    match human_chromosome_index(chromosome) {
        Some(i) => {
            if ref_genome_is_v37() {
                CENTROMERES_37[i]
            } else {
                CENTROMERES_38[i]
            }
        }
        None => 0,
    }
}

/// SvUtilities.getChromosomalArm
pub fn get_chromosomal_arm(chromosome: &str, position: i32) -> Arm {
    if position < get_centromere(chromosome) { Arm::P } else { Arm::Q }
}

/// SvUtilities.getChromosomalArmLength
pub fn get_chromosomal_arm_length(chromosome: &str, arm: Arm) -> i32 {
    // HumanChromosome.fromString throws for other contigs
    let i = human_chromosome_index(chromosome).expect("human chromosome");
    let centromere = if ref_genome_is_v37() { CENTROMERES_37[i] } else { CENTROMERES_38[i] };
    if arm == Arm::P {
        return centromere;
    }
    get_chromosome_length(chromosome) - centromere
}

/// SvUtilities.addSvToChrBreakendMap: each breakend into its chromosome's list by ascending position
pub fn add_sv_to_chr_breakend_map(var: &Sv, map: &mut JavaHashMap<String, Vec<Be>>) {
    for be in [SE_START, SE_END] {
        if be == SE_END && var.is_sgl_breakend() {
            continue;
        }
        let breakend = var.breakend_se(be).unwrap();
        let position = breakend.position();
        let chr = breakend.chromosome().to_string();
        if !map.contains_key(&chr) {
            map.insert(chr.clone(), Vec::new());
        }
        let list = map.get_mut(&chr).unwrap();
        let mut index = 0;
        while index < list.len() {
            let other = &list[index];
            if position < other.position() {
                break;
            }
            // an inferred placed at another SV's location faces its partner breakend
            if position == other.position() && breakend.orientation() != other.orientation() {
                if var.is_inferred_sgl() && breakend.reverse_orient() {
                    break;
                } else if other.sv().is_inferred_sgl() && other.forward_orient() {
                    break;
                }
            }
            index += 1;
        }
        list.insert(index, breakend);
    }
}

/// SvUtilities.findCentromereBreakendIndex
pub fn find_centromere_breakend_index(list: Option<&Vec<Be>>, arm: Arm) -> i32 {
    let Some(list) = list else { return -1 };
    if list.is_empty() {
        return -1;
    }
    if arm == Arm::P {
        let mut i = 0;
        while i < list.len() {
            if list[i].arm() == Arm::Q {
                break;
            }
            i += 1;
        }
        i as i32 - 1
    } else {
        let mut i = list.len() as i32 - 1;
        while i >= 0 {
            if list[i as usize].arm() == Arm::P {
                break;
            }
            i -= 1;
        }
        if i + 1 >= list.len() as i32 { -1 } else { i + 1 }
    }
}

/// SvUtilities.getSvTypesStr: "BND=1_DEL=2" in type order
pub fn get_sv_types_str(counts: &[i32; 7]) -> String {
    let mut s = String::new();
    for (i, c) in counts.iter().enumerate() {
        if *c == 0 {
            continue;
        }
        let item = format!("{}={}", hmf_common::sv::StructuralVariantType::VALUES[i], c);
        append_str(&mut s, &item, '_');
    }
    s
}

pub fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}

/// SvUtilities.isWithin
pub fn is_within(variant: &SvVarData, chromosome: &str, position: i32) -> bool {
    if variant.chromosome(true) != chromosome || variant.chromosome(false) != chromosome {
        return false;
    }
    position_within(position, variant.position(true), variant.position(false))
}

/// SvUtilities.isOverlapping
pub fn is_overlapping(v1: &SvVarData, v2: &SvVarData) -> bool {
    is_within(v2, v1.chromosome(true), v1.position(true))
        || is_within(v2, v1.chromosome(false), v1.position(false))
        || is_within(v1, v2.chromosome(true), v2.position(true))
        || is_within(v1, v2.chromosome(false), v2.position(false))
}

/// SvUtilities.getProximity
pub fn get_proximity(var1: &SvVarData, var2: &SvVarData) -> i32 {
    let mut min_distance = -1;
    for se1 in [SE_START, SE_END] {
        let Some(be1) = var1.breakend_se(se1) else { continue };
        for se2 in [SE_START, SE_END] {
            let Some(be2) = var2.breakend_se(se2) else { continue };
            if be1.chromosome() == be2.chromosome() {
                let d = (be1.position() - be2.position()).abs();
                if min_distance == -1 || d < min_distance {
                    min_distance = d;
                }
            }
        }
    }
    min_distance
}

pub fn calc_consistency_list(svs: &[Sv]) -> i32 {
    svs.iter().map(|v| calc_consistency_sv(v)).sum()
}
pub fn calc_consistency_sv(var: &SvVarData) -> i32 {
    let mut c = 0;
    for be in [SE_START, SE_END] {
        if be == SE_END && var.is_sgl_breakend() {
            continue;
        }
        c += calc_consistency_be(&var.breakend_se(be).unwrap());
    }
    c
}
pub fn calc_consistency_be(breakend: &SvBreakend) -> i32 {
    (if breakend.arm() == Arm::P { 1 } else { -1 }) * breakend.orientation() as i32
}

pub const MAX_COPY_NUM_DIFF: f64 = 0.5;
pub const MAX_COPY_NUM_DIFF_PERC: f64 = 0.15;

pub fn copy_numbers_equal(cn1: f64, cn2: f64) -> bool {
    copy_numbers_equal_with(cn1, cn2, MAX_COPY_NUM_DIFF, MAX_COPY_NUM_DIFF_PERC)
}
pub fn copy_numbers_equal_with(cn1: f64, cn2: f64, max_diff: f64, max_diff_perc: f64) -> bool {
    let diff = (cn2 - cn1).abs();
    let perc = diff / jmax(cn1.abs(), cn2.abs());
    !(diff > max_diff && perc > max_diff_perc)
}

/// SvUtilities.formatJcn (for logging)
pub fn format_jcn(jcn: f64) -> String {
    if jcn > 10.0 {
        hmf_common::fmt::java_fixed(jcn, 0)
    } else if jcn < 0.5 {
        hmf_common::fmt::java_fixed(jcn, 2)
    } else {
        hmf_common::fmt::java_fixed(jcn, 1)
    }
}

pub fn make_chr_arm_str_sv(var: &SvVarData, use_start: bool) -> String {
    make_chr_arm_str(var.chromosome(use_start), var.arm(use_start))
}

/// SvUtilities.loadConfigFile: regions from resource lines (chromosome, start, end), skipping header lines
pub fn load_config_file(lines: &[String], v37: Option<bool>, delim: char) -> Option<Vec<hmf_common::regions::ChrBaseRegion>> {
    let mut regions = Vec::new();
    for line in lines {
        if line.contains("Chromosome") {
            continue;
        }
        let items: Vec<&str> = line.split(delim).collect();
        if items.len() < 3 {
            hmf_common::error!("invalid resource file: {}", line);
            return None;
        }
        let chromosome = match v37 {
            Some(v) => hmf_common::regions::versioned_chromosome(!v, items[0]),
            None => items[0].to_string(),
        };
        regions.push(hmf_common::regions::ChrBaseRegion::new(&chromosome, items[1].parse().ok()?, items[2].parse().ok()?));
    }
    Some(regions)
}

/// Java's `(int)` cast of a double: truncation towards zero, NaN to 0, saturating
pub fn d2i(d: f64) -> i32 {
    d as i32
}

/// String.format("%.Nf")
pub fn ff(d: f64, places: usize) -> String {
    hmf_common::fmt::java_fixed(d, places)
}

/// String.valueOf(double)
pub fn dstr(d: f64) -> String {
    hmf_common::fmt::java_double_to_string(d)
}
