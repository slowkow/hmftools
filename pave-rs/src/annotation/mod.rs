//! The annotation resources: PON, gnomAD frequencies, mappability, ClinVar and the blacklists.
//!
//! Java loads each resource a chromosome at a time and keeps all of a chromosome's entries in memory
//! (the PON, gnomAD and mappability files are 0.8-2.7 GB of text). The port reads the input VCF first,
//! so it knows the positions its variants need, and keeps only the entries at those positions; each
//! lookup then gives what Java's would (the PON, gnomAD and ClinVar lookups match exact positions,
//! mappability takes the entry at or before a position).

pub mod blacklist;
pub mod clinvar;
pub mod gnomad;
pub mod mappability;
pub mod pon;
pub mod reportability;

use std::collections::HashMap;

/// The positions the variants need, by chromosome, sorted and without repeats
#[derive(Default, Clone)]
pub struct NeededPositions {
    pub by_chr: HashMap<String, Vec<i32>>,
}

impl NeededPositions {
    pub fn get(&self, chromosome: &str) -> Option<&[i32]> {
        self.by_chr.get(chromosome).map(|v| v.as_slice())
    }

    pub fn finish(&mut self) {
        for v in self.by_chr.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
    }

    #[allow(dead_code)]
    pub fn contains(&self, chromosome: &str, position: i32) -> bool {
        self.get(chromosome).is_some_and(|v| v.binary_search(&position).is_ok())
    }
}

/// Integer.parseInt of a field as bytes (plain digits parsed directly, anything else as Java parses it)
#[inline]
pub fn parse_int(field: &[u8]) -> Option<i32> {
    if !field.is_empty() && field.len() <= 9 && field.iter().all(|b| b.is_ascii_digit()) {
        return Some(field.iter().fold(0i32, |n, b| n * 10 + (b - b'0') as i32));
    }
    hmf_common::java_text::parse_int(std::str::from_utf8(field).ok()?)
}

/// Double.parseDouble of a field as bytes
pub fn parse_double(field: &[u8]) -> Option<f64> {
    hmf_common::java_text::parse_double(std::str::from_utf8(field).ok()?)
}

/// The fields of a line split at a delimiter (String.split(delim, -1): every field kept)
pub fn split_fields(line: &[u8], delim: u8) -> Vec<&[u8]> {
    line.split(|b| *b == delim).collect()
}

/// CommonInfo.putAttribute(key, value) without allowing overwrites, as PAVE puts CLNSIG, CLNSIGCONF, the
/// blacklist flags, HOTSPOT and REPORTED.
///
/// PA9: htsjdk throws if the record has the key already (an input annotated by PAVE before, or a SAGE
/// record with the HOTSPOT flag under -set_reportable), which ends the run. Corrected: the value is
/// replaced, as PAVE does for PAVE_TI and IMPACT.
pub fn put_once(rec: &mut hmf_common::vcf::VcfRecord, key: &str, value: hmf_common::vcf::InfoValue) -> anyhow::Result<()> {
    if rec.has_attr(key) && crate::compat::java("PA9") {
        anyhow::bail!("IllegalStateException: Attempting to overwrite key->value binding: key = {key}");
    }
    rec.put_attr(key, value);
    Ok(())
}

pub fn text(field: &[u8]) -> String {
    String::from_utf8_lossy(field).into_owned()
}
