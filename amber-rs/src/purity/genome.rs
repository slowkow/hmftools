//! Reference-genome facts the tumor-only steps take from hmf-common: HumanChromosome's order, the centromere
//! positions and centromere gaps (RefGenomeCoordinates, from the resources refgenome/centromeres.<v>.tsv and
//! refgenome/centromere_gaps.<v>.tsv), and the HLA, IG and TR regions (ImmuneRegions).
//!
//! Shared by the purity analysis and LohIdentifier (crate::loh); it could move to the crate root or to
//! hmf-common-rs when the copies are merged.

use hmf_common::regions::ChrBaseRegion;

/// HumanChromosome.fromString's ordinal (0 for chromosome 1 ... 21 for 22, 22 for X, 23 for Y), or None for
/// a contig that is not a human chromosome (where Java throws IllegalArgumentException). Java strips a
/// leading "chr" in any case.
pub fn chromosome_ordinal(chromosome: &str) -> Option<usize> {
    let bare = if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome };
    match bare {
        "X" => Some(22),
        "Y" => Some(23),
        _ => {
            // HumanChromosome.valueOf("_" + name): the enum names are _1 to _22, no leading zeros or signs
            if bare.is_empty() || bare.len() > 2 || !bare.bytes().all(|b| b.is_ascii_digit()) || bare.starts_with('0') {
                return None;
            }
            let n: usize = bare.parse().ok()?;
            if (1..=22).contains(&n) { Some(n - 1) } else { None }
        }
    }
}

/// HumanChromosome.isAllosome: X and Y
pub fn is_allosome(chromosome: &str) -> bool {
    matches!(chromosome_ordinal(chromosome), Some(22) | Some(23))
}

const CENTROMERES_37: [i32; 24] = [
    123035434, 93826171, 92004854, 51160117, 47905641, 60330166, 59554331, 45338887, 48867679, 40754935, 53144205, 36356694, 17500000,
    17500000, 18500000, 36835801, 23763006, 16960898, 26181782, 27869569, 12788129, 14500000, 60132012, 11604553,
];

const CENTROMERES_38: [i32; 24] = [
    123605523, 93139351, 92214016, 50726026, 48272854, 59191911, 59498944, 44955505, 44377363, 40640102, 52751711, 35977330, 17025624,
    17086762, 18362627, 37295920, 24849830, 18161053, 25844927, 28237290, 11890184, 14004553, 60509061, 10430492,
];

/// RefGenomeCoordinates.centromere: 0 for a contig that HumanChromosome.contains rejects
pub fn centromere(v38: bool, chromosome: &str) -> i32 {
    // HumanChromosome.contains strips a lower-case "chr" only; fromString is used for the lookup
    let bare = chromosome.strip_prefix("chr").unwrap_or(chromosome);
    let contained = match bare {
        "X" | "Y" => true,
        _ => bare.bytes().all(|b| b.is_ascii_digit()) && !bare.is_empty() && bare.parse::<i64>().is_ok_and(|n| (1..=22).contains(&n)),
    };
    if !contained {
        return 0;
    }
    match chromosome_ordinal(chromosome) {
        Some(i) => {
            if v38 {
                CENTROMERES_38[i]
            } else {
                CENTROMERES_37[i]
            }
        }
        // "chr01" is contained (Integer.parseInt) but fromString fails: Java throws
        None => panic!("No enum constant com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._{bare}"),
    }
}

/// refgenome/centromere_gaps.37.tsv: (chromosome ordinal, start, end) in file order
const CENTROMERE_GAPS_37: [(usize, i32, i32); 31] = [
    (0, 121535434, 142535434),
    (1, 92326171, 95326171),
    (2, 90504854, 93504854),
    (3, 49660117, 52660117),
    (4, 46405641, 49405641),
    (5, 58830166, 61830166),
    (6, 58054331, 61054331),
    (7, 43838887, 46838887),
    (8, 47367679, 65367679),
    (9, 39254935, 42254935),
    (10, 51644205, 54644205),
    (11, 34856694, 37856694),
    (12, 10000, 19000000),
    (13, 10000, 19000000),
    (14, 10000, 20000000),
    (15, 35285801, 46385801),
    (16, 22263006, 25263006),
    (17, 15460898, 18460898),
    (18, 24681782, 27681782),
    (19, 26369569, 29369569),
    (20, 0, 9411193),
    (20, 9595548, 9645548),
    (20, 9775437, 9825437),
    (20, 10034920, 10084920),
    (20, 10215976, 10365976),
    (20, 10647896, 10697896),
    (20, 11188129, 14338129),
    (21, 10000, 16000000),
    (22, 58632012, 61632012),
    (23, 10104553, 13104553),
    (23, 28819361, 58819361),
];

/// refgenome/centromere_gaps.38.tsv: (chromosome ordinal, start, end) in file order
const CENTROMERE_GAPS_38: [(usize, i32, i32); 33] = [
    (0, 125184587, 143184587),
    (8, 45518558, 60518558),
    (23, 26673214, 56673214),
    (12, 10000, 16000000),
    (13, 10000, 16000000),
    (14, 10000, 17000000),
    (15, 38280682, 46280682),
    (20, 10000, 5010000),
    (21, 10000, 10510000),
    (0, 123105523, 124105523),
    (1, 92839351, 93639351),
    (2, 91914016, 92714016),
    (3, 50426026, 51226026),
    (4, 47972854, 48772854),
    (5, 58891911, 59691911),
    (6, 59198944, 59998944),
    (7, 44655505, 45455505),
    (8, 44077363, 44877363),
    (9, 40340102, 41140102),
    (10, 52451711, 53251711),
    (11, 35677330, 36477330),
    (12, 16725624, 17525624),
    (13, 16786762, 17586762),
    (14, 18062627, 18862627),
    (15, 36995920, 37795920),
    (16, 24549830, 25349830),
    (17, 17861053, 18661053),
    (18, 25544927, 26344927),
    (19, 27937290, 28737290),
    (20, 11590184, 12390184),
    (21, 13704553, 14504553),
    (22, 60209061, 61009061),
    (23, 10130492, 10930492),
];

/// RefGenomeCoordinates.centromereGaps(version).get(chromosome): the gaps of one chromosome as (start, end),
/// in file order (empty where Java's map has no entry)
pub fn centromere_gaps(v38: bool, chromosome: &str) -> Vec<(i32, i32)> {
    let Some(ordinal) = chromosome_ordinal(chromosome) else { return Vec::new() };
    let gaps: &[(usize, i32, i32)] = if v38 { &CENTROMERE_GAPS_38 } else { &CENTROMERE_GAPS_37 };
    gaps.iter().filter(|g| g.0 == ordinal).map(|g| (g.1, g.2)).collect()
}

const TRANSCRIPT_FLANK: i32 = 5_000;

/// ImmuneRegions.getHlaRegions: HLA-A, HLA-B, HLA-C with 5kb flanks
pub fn hla_regions(v38: bool) -> Vec<ChrBaseRegion> {
    if v38 {
        vec![
            ChrBaseRegion::new("chr6", 29942532 - TRANSCRIPT_FLANK, 29945870 + TRANSCRIPT_FLANK),
            ChrBaseRegion::new("chr6", 31353872 - TRANSCRIPT_FLANK, 31367067 + TRANSCRIPT_FLANK),
            ChrBaseRegion::new("chr6", 31268749 - TRANSCRIPT_FLANK, 31272092 + TRANSCRIPT_FLANK),
        ]
    } else {
        vec![
            ChrBaseRegion::new("6", 29909037 - TRANSCRIPT_FLANK, 29913661 + TRANSCRIPT_FLANK),
            ChrBaseRegion::new("6", 31321649 - TRANSCRIPT_FLANK, 31324964 + TRANSCRIPT_FLANK),
            ChrBaseRegion::new("6", 31236526 - TRANSCRIPT_FLANK, 31239863 + TRANSCRIPT_FLANK),
        ]
    }
}

/// ImmuneRegions.getIgRegions: IGK, IGH, IGL
pub fn ig_regions(v38: bool) -> Vec<ChrBaseRegion> {
    if v38 {
        vec![
            ChrBaseRegion::new("chr2", 88857161, 90315836),
            ChrBaseRegion::new("chr14", 105586437, 106879844),
            ChrBaseRegion::new("chr22", 22026076, 22922913),
        ]
    } else {
        vec![
            ChrBaseRegion::new("2", 89156674, 90538397),
            ChrBaseRegion::new("14", 106032614, 107288051),
            ChrBaseRegion::new("22", 22380474, 23265085),
        ]
    }
}

/// ImmuneRegions.getTrRegions: TRG, TRB, TRA/TRD
pub fn tr_regions(v38: bool) -> Vec<ChrBaseRegion> {
    if v38 {
        vec![
            ChrBaseRegion::new("chr7", 38239580, 38367882),
            ChrBaseRegion::new("chr7", 142299177, 142812869),
            ChrBaseRegion::new("chr14", 21622293, 22552156),
        ]
    } else {
        vec![
            ChrBaseRegion::new("7", 38279181, 38407483),
            ChrBaseRegion::new("7", 141999017, 142510554),
            ChrBaseRegion::new("14", 22090446, 23021099),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chromosomes() {
        assert_eq!(chromosome_ordinal("chr1"), Some(0));
        assert_eq!(chromosome_ordinal("1"), Some(0));
        assert_eq!(chromosome_ordinal("CHR22"), Some(21));
        assert_eq!(chromosome_ordinal("chrX"), Some(22));
        assert_eq!(chromosome_ordinal("Y"), Some(23));
        assert_eq!(chromosome_ordinal("chrM"), None);
        assert_eq!(chromosome_ordinal("23"), None);
        assert_eq!(chromosome_ordinal("01"), None);
        assert!(is_allosome("chrX") && is_allosome("Y") && !is_allosome("chr7"));
        assert_eq!(centromere(true, "chr1"), 123605523);
        assert_eq!(centromere(false, "Y"), 11604553);
        assert_eq!(centromere(true, "chrM"), 0);
    }

    #[test]
    fn gaps() {
        assert_eq!(centromere_gaps(true, "chr1"), vec![(125184587, 143184587), (123105523, 124105523)]);
        assert_eq!(centromere_gaps(true, "chrY"), vec![(26673214, 56673214), (10130492, 10930492)]);
        assert_eq!(centromere_gaps(false, "21").len(), 7);
        assert_eq!(centromere_gaps(false, "Y"), vec![(10104553, 13104553), (28819361, 58819361)]);
        assert!(centromere_gaps(true, "chrM").is_empty());
    }
}
