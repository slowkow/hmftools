//! HgvsCoding: the HGVS coding (c.) or non-coding (n.) string of a variant in a transcript.
//!
//! Conventions: coding bases on the negative strand are reversed; the start codon is position 1;
//! -1 is the base before the start codon and *1 the first base after the stop codon; the most 3'
//! position is used for duplications and deletions in homology; intronic positions add the distance
//! into the intron (+ after an exon, - before one); an insertion is a duplication when the
//! microhomology matches the inserted bases.

use hmf_common::codon::reverse_complement_bases;

use super::effects::{CodingType, RegionType};
use super::{CodingContext, NEG_STRAND, POS_STRAND, SE_END, SE_START};
use crate::variant::VariantData;

/// HGVS_UNKNOWN: a misclassified or malformed string
pub const HGVS_UNKNOWN: &str = "unknown";

pub fn set(variant: &VariantData, cc: &mut CodingContext) {
    cc.hgvs = generate(variant, cc);
}

pub fn generate(variant: &VariantData, cc: &CodingContext) -> String {
    if cc.region_type == RegionType::Upstream && cc.coding_type != CodingType::Enhancer {
        return String::new(); // undefined by the standard
    }
    let mut sb = String::with_capacity(32);
    sb.push_str(if cc.coding_type == CodingType::NonCoding { "n." } else { "c." });
    if variant.is_base_change() || variant.is_mixed() {
        form_point_mutation(variant, cc, &mut sb);
    } else if variant.is_deletion() {
        form_deletion(variant, cc, &mut sb);
    } else if variant.is_insert() {
        form_insertion(variant, cc, &mut sb);
    }
    sb
}

fn add_utr(cc: &CodingContext, sb: &mut String) {
    if cc.coding_ends_on_exon_boundary {
        return;
    }
    if cc.coding_type == CodingType::Utr5p || cc.coding_type == CodingType::Enhancer {
        sb.push('-');
    } else if cc.coding_type == CodingType::Utr3p {
        sb.push('*');
    }
}

fn add_coding_base(cc: &CodingContext, mut coding_base: i32, sb: &mut String, is_second: bool) {
    if is_second {
        sb.push('_');
    }
    if is_second && cc.spans_coding_end {
        return;
    }
    add_utr(cc, sb);
    if cc.coding_base == 0 {
        sb.push('-');
        coding_base = 1;
    }
    sb.push_str(&coding_base.to_string());
}

fn add_intronic_position(position: i32, spans_stop_codon: bool, sb: &mut String) {
    if spans_stop_codon {
        // MNVs or DELs which delete the stop codon
        sb.push('*');
    } else if position > 0 {
        sb.push('+');
    } else {
        sb.push('-');
    }
    sb.push_str(&position.abs().to_string());
}

fn add_bases(bases: &str, cc: &CodingContext, sb: &mut String) {
    if cc.strand == POS_STRAND {
        sb.push_str(bases);
    } else {
        sb.push_str(&String::from_utf8(reverse_complement_bases(bases.as_bytes())).unwrap());
    }
}

/// eg c.191delG, c.203_204delCG, c.60+5601_60+5610delTTTTTTTTTT, c.123-7190_123-7175del..., a splice
/// region overlap c.10_12+3delAGTCCC or c.7_9-2delGATCC, across the coding start c-2_4delATGC
fn form_deletion(variant: &VariantData, cc: &CodingContext, sb: &mut String) {
    let del_length = variant.base_diff().abs();
    let coding_bases_deleted = cc.deleted_coding_bases;
    let coding_base = cc.coding_base;
    let mut nearest_exon = cc.nearest_exon_distance;

    // the variant spanning into the intron
    let spans_coding_boundary = cc.region_type == RegionType::Exonic && del_length > coding_bases_deleted;
    let spans_coding_start = cc.spans_coding_start && nearest_exon < -1; // cannot be just the 1st 5'UTR base
    let spans_splice = !cc.spans_coding_start && spans_coding_boundary && nearest_exon != 0;
    let spans_upstream_splice = spans_splice && nearest_exon < 0;
    let spans_downstream_splice = spans_splice && !spans_upstream_splice;

    // if the deleted bases don't cross a splice region but the range was trimmed to an exon or coding
    // boundary, don't move past the first (non-deleted) ref base
    let coding_base_range_capped = (cc.strand == POS_STRAND && cc.coding_position_range[SE_START] > variant.position)
        || (cc.strand == NEG_STRAND && cc.coding_position_range[SE_END] < variant.end_position);

    let coding_base_start;
    let coding_base_end;

    if cc.region_type == RegionType::Exonic {
        if spans_coding_start {
            coding_base_start = nearest_exon + 1;
            nearest_exon = 0; // not intronic
            coding_base_end = coding_bases_deleted;
        } else if cc.coding_type == CodingType::Utr5p {
            // the furthest from the start of coding, then back one from the ref position
            coding_base_start = if !spans_upstream_splice { coding_base - 1 } else { coding_base };
            coding_base_end = coding_base_start - coding_bases_deleted + 1;
        } else {
            // move past the ref base
            coding_base_start = if !spans_upstream_splice && !coding_base_range_capped { coding_base + 1 } else { coding_base };
            coding_base_end = if del_length > 1 { coding_base_start + coding_bases_deleted - 1 } else { coding_base_start };
        }
    } else {
        coding_base_start = coding_base;
        coding_base_end = coding_base;
    }

    nearest_exon += 1;

    add_coding_base(cc, coding_base_start, sb, false);

    if cc.region_type == RegionType::Intronic || spans_upstream_splice {
        add_intronic_position(nearest_exon, false, sb);
    }

    if del_length > 1 {
        add_coding_base(cc, coding_base_end, sb, true);
        if cc.region_type == RegionType::Intronic || spans_downstream_splice {
            if spans_downstream_splice {
                nearest_exon = cc.nearest_exon_distance - 1;
            } else {
                nearest_exon += del_length - 1;
            }
            add_intronic_position(nearest_exon, cc.spans_coding_end, sb);
        }
    }

    sb.push_str("del");
    add_bases(&variant.ref_[1..], cc, sb);
}

pub fn is_duplication(variant: &VariantData) -> bool {
    variant.is_insert() && variant.microhomology() == &variant.alt[1..]
}

/// eg c.1033_1034insA, c.15+1619_15+1620insTTTGTT, c.-23-304_-23-303insA
fn form_insertion(variant: &VariantData, cc: &CodingContext, sb: &mut String) {
    if is_duplication(variant) {
        form_duplication(variant, cc, sb);
        return;
    }
    let inserted = &variant.alt[1..];
    let mut coding_base = cc.coding_base;
    let mut nearest_exon = cc.nearest_exon_distance;

    add_coding_base(cc, coding_base, sb, false);
    if cc.region_type == RegionType::Intronic {
        add_intronic_position(nearest_exon, false, sb);
    }
    if cc.region_type == RegionType::Exonic {
        if cc.coding_type == CodingType::Utr5p {
            coding_base -= 1;
        } else {
            coding_base += 1;
        }
    }
    add_coding_base(cc, coding_base, sb, true);
    if cc.region_type == RegionType::Intronic {
        nearest_exon += 1;
        add_intronic_position(nearest_exon, false, sb);
    }
    sb.push_str("ins");
    add_bases(inserted, cc, sb);
}

/// eg c.377+9769_377+9772dupAAAT, c.1440-82dupT (a single location for a single base)
fn form_duplication(variant: &VariantData, cc: &CodingContext, sb: &mut String) {
    let mut coding_base = cc.coding_base;

    if cc.region_type == RegionType::Exonic && cc.strand == POS_STRAND {
        // move forward to the start of the duplicated section; the negative strand's coding base is on the
        // upper side of the INDEL already
        if variant.is_realigned_variant() {
            // walk back to the start of the repeat
            let mc_length = variant.microhomology().len() as i32;
            coding_base -= mc_length - 1;
        } else {
            coding_base += 1;
        }
    }

    // intronic bases with a +ve nearest exon are shifted by 1 to the start of the duplication likewise
    let nearest_exon = cc.nearest_exon_distance;
    let intron_base_adjustment = if cc.strand == POS_STRAND { 1 } else { 0 };

    let inserted = &variant.alt[1..];
    let insert_length = inserted.len() as i32;

    if insert_length == 1 {
        add_coding_base(cc, coding_base, sb, false);
        if cc.region_type == RegionType::Intronic {
            add_intronic_position(nearest_exon + intron_base_adjustment, false, sb);
        }
    } else {
        let mut base_shift = insert_length - 1;
        if cc.coding_type == CodingType::Utr5p {
            base_shift = -base_shift;
        }
        // on the negative strand the positions are flipped around
        let (coding_base_start, coding_base_end, mut intron_base_start, mut intron_base_end);
        if cc.strand == POS_STRAND {
            coding_base_start = coding_base;
            coding_base_end = if cc.region_type == RegionType::Exonic { coding_base + base_shift } else { coding_base };
            intron_base_start = nearest_exon;
            intron_base_end = if cc.region_type == RegionType::Intronic { nearest_exon + base_shift } else { nearest_exon };
        } else {
            coding_base_end = coding_base;
            coding_base_start = if cc.region_type == RegionType::Exonic { coding_base - base_shift } else { coding_base };
            intron_base_end = nearest_exon;
            intron_base_start = if cc.region_type == RegionType::Intronic { nearest_exon - base_shift } else { nearest_exon };
        }
        intron_base_start += intron_base_adjustment;
        intron_base_end += intron_base_adjustment;

        add_coding_base(cc, coding_base_start, sb, false);
        if cc.region_type == RegionType::Intronic {
            add_intronic_position(intron_base_start, false, sb);
        }
        add_coding_base(cc, coding_base_end, sb, true);
        if cc.region_type == RegionType::Intronic {
            add_intronic_position(intron_base_end, false, sb);
        }
    }

    sb.push_str("dup");
    add_bases(inserted, cc, sb);
}

/// eg c.-14G>C, c.76A>C, c.88+1G>T, c.89-2A>C; MNVs c.38-20810_38-20808delTTGinsACT, c.132_133delCCinsTT,
/// c.*907_*908delCCinsGT
fn form_point_mutation(variant: &VariantData, cc: &CodingContext, sb: &mut String) {
    let coding_base = cc.coding_base;
    let mut nearest_exon = cc.nearest_exon_distance;
    let var_length = variant.ref_.len().max(variant.alt.len()) as i32;

    if var_length == 1 {
        add_coding_base(cc, coding_base, sb, false);
        if cc.region_type == RegionType::Intronic {
            add_intronic_position(nearest_exon, false, sb);
        }
        add_bases(&variant.ref_, cc, sb);
        sb.push('>');
        add_bases(&variant.alt, cc, sb);
    } else {
        let spans_splice = !cc.spans_coding_start && cc.region_type == RegionType::Exonic && nearest_exon != 0;
        let spans_upstream_splice = spans_splice && nearest_exon < 0;
        let spans_downstream_splice = spans_splice && !spans_upstream_splice;

        let mut coding_base_start = coding_base;
        let mut coding_base_len = if cc.coding_type == CodingType::Enhancer {
            var_length
        } else {
            cc.coding_position_range[SE_END] - cc.coding_position_range[SE_START] + 1
        };
        if variant.is_mixed() {
            // the last coding base matches the ref for these variants
            coding_base_len -= 1;
        }
        let mut coding_base_end = if cc.region_type == RegionType::Exonic { coding_base + coding_base_len - 1 } else { coding_base };

        if (cc.coding_type == CodingType::Utr5p && cc.region_type == RegionType::Exonic) || cc.coding_type == CodingType::Enhancer {
            // the highest exonic base (furthest from the coding start) was recorded; the coding bases are
            // negated but drop from the highest to lowest, eg c.-83_-80delGGATinsCTGG
            coding_base_start = coding_base;
            coding_base_end = coding_base - coding_base_len + 1;
        } else if cc.spans_coding_start {
            coding_base_start = nearest_exon;
            nearest_exon = 0; // not intronic
            coding_base_end = coding_base_len;
        }

        add_coding_base(cc, coding_base_start, sb, false);
        if cc.region_type == RegionType::Intronic || spans_upstream_splice {
            add_intronic_position(nearest_exon, false, sb);
        }

        if !variant.is_insert() {
            // for mixed INDELs
            add_coding_base(cc, coding_base_end, sb, true);
            if cc.region_type == RegionType::Intronic || spans_downstream_splice {
                if spans_downstream_splice {
                    nearest_exon = cc.nearest_exon_distance;
                } else {
                    nearest_exon += var_length - 1;
                }
                add_intronic_position(nearest_exon, cc.spans_coding_end, sb);
            }
        }

        sb.push_str("del");
        add_bases(&variant.ref_, cc, sb);
        sb.push_str("ins");
        add_bases(&variant.alt, cc, sb);
    }
}
