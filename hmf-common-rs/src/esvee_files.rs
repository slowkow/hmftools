//! hmf-common `sv.EsveeDiscordantStats` (added 2026-11-27 for qsee): ESVEE prep's discordant fragment
//! statistics (`.esvee.disc_stats.tsv`, formerly `.esvee.prep.disc_stats.tsv`).

use anyhow::{anyhow, Context, Result};

use crate::java_text::{parse_long, read_lines, split};

/// DiscordantFragType, in enum order
pub const DISCORDANT_FRAG_TYPES: [&str; 11] = [
    "Translocation",
    "InvLt1K",
    "Inv1To5K",
    "Inv5To100K",
    "InvGt100K",
    "Del1To5K",
    "Del5To100K",
    "DelGt100K",
    "Dup1To5K",
    "Dup5To100K",
    "DupGt100K",
];

/// EsveeDiscordantStats
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EsveeDiscordantStats {
    pub total_reads: i64,
    pub prep_reads: i64,
    /// counts in DiscordantFragType order
    pub type_counts: [i64; 11],
}

/// EsveeDiscordantStats.generateFilename: basePath + File.separator + sample + ".esvee.disc_stats.tsv"
pub fn esvee_disc_stats_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.esvee.disc_stats.tsv")
}

/// EsveeDiscordantStats.read
pub fn read_esvee_disc_stats(path: &str) -> Result<EsveeDiscordantStats> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    esvee_disc_stats_from_lines(&lines).with_context(|| format!("reading {path}"))
}

/// EsveeDiscordantStats.read on lines already read: the second line, split at tabs (trailing empty
/// fields dropped): total reads, prep reads, then the counts by type
pub fn esvee_disc_stats_from_lines(lines: &[String]) -> Result<EsveeDiscordantStats> {
    let line = lines.get(1).ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index 1 out of bounds for length {}", lines.len()))?;
    let v = split(line, '\t');
    let at = |i: usize| -> Result<i64> {
        let s = v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))?;
        parse_long(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
    };
    let total_reads = at(0)?;
    let prep_reads = at(1)?;
    let mut type_counts = [0i64; 11];
    for (i, c) in type_counts.iter_mut().enumerate() {
        *c = at(i + 2)?;
    }
    Ok(EsveeDiscordantStats { total_reads, prep_reads, type_counts })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disc_stats_lines() {
        let lines: Vec<String> = ["TotalReads\tPrepReads\tTranslocation", "100\t10\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11"].iter().map(|s| s.to_string()).collect();
        let s = esvee_disc_stats_from_lines(&lines).unwrap();
        assert_eq!((s.total_reads, s.prep_reads, s.type_counts[0], s.type_counts[10]), (100, 10, 1, 11));
        assert!(esvee_disc_stats_from_lines(&lines[..1]).is_err());
        assert!(esvee_disc_stats_from_lines(&[lines[0].clone(), "1\t2\t3".to_string()]).is_err());
    }
}
