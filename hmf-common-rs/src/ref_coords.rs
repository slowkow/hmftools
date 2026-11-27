//! hmf-common RefGenomeCoordinates: the lengths of the human chromosomes in GRCh37 and GRCh38
//! (resources refgenome/lengths.37.tsv and lengths.38.tsv).

/// 1-22, X, Y
const NAMES: [&str; 24] = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"];

const LENGTHS_37: [i32; 24] = [
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
];

const LENGTHS_38: [i32; 24] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
];

/// The length of a human chromosome given without a "chr" prefix ("1".."22", "X", "Y"; as
/// HumanChromosome's enum names, case-sensitive), or None
pub fn chromosome_length(chromosome: &str, v37: bool) -> Option<i32> {
    let i = NAMES.iter().position(|n| *n == chromosome)?;
    Some(if v37 { LENGTHS_37[i] } else { LENGTHS_38[i] })
}

/// RefGenomeCoordinates' centromeres (resources refgenome/centromeres.37.tsv and centromeres.38.tsv), in the
/// order of NAMES. (Added for the LINX visualiser.)
const CENTROMERES_37: [i32; 24] = [
    123035434, 93826171, 92004854, 51160117, 47905641, 60330166, 59554331, 45338887, 48867679, 40754935, 53144205, 36356694,
    17500000, 17500000, 18500000, 36835801, 23763006, 16960898, 26181782, 27869569, 12788129, 14500000, 60132012, 11604553,
];

const CENTROMERES_38: [i32; 24] = [
    123605523, 93139351, 92214016, 50726026, 48272854, 59191911, 59498944, 44955505, 44377363, 40640102, 52751711, 35977330,
    17025624, 17086762, 18362627, 37295920, 24849830, 18161053, 25844927, 28237290, 11890184, 14004553, 60509061, 10430492,
];

/// The centromere position of a human chromosome given without a "chr" prefix, or None.
pub fn centromere(chromosome: &str, v37: bool) -> Option<i32> {
    let i = NAMES.iter().position(|n| *n == chromosome)?;
    Some(if v37 { CENTROMERES_37[i] } else { CENTROMERES_38[i] })
}

/// RefGenomeCoordinates.readCentromereGaps: the lines of refgenome/centromere_gaps.<version>.tsv
/// ("hs<chromosome>", start, end; tab-separated), as BufferedReader.lines() gives them.
pub fn centromere_gap_lines(v37: bool) -> Vec<&'static str> {
    let text = if v37 { include_str!("../resources/refgenome/centromere_gaps.37.tsv") } else { include_str!("../resources/refgenome/centromere_gaps.38.tsv") };
    text.lines().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lengths() {
        assert_eq!(chromosome_length("1", false), Some(248956422));
        assert_eq!(chromosome_length("X", true), Some(155270560));
        assert_eq!(chromosome_length("Y", false), Some(57227415));
        assert_eq!(chromosome_length("M", false), None);
        assert_eq!(chromosome_length("x", false), None);
    }

    #[test]
    fn centromeres() {
        assert_eq!(centromere("1", false), Some(123605523));
        assert_eq!(centromere("Y", true), Some(11604553));
        assert_eq!(centromere("chr1", false), None);
        assert_eq!(centromere_gap_lines(false)[0], "hs1\t125184587\t143184587");
        // the GRCh38 file ends without a newline
        assert_eq!(centromere_gap_lines(false).len(), 33);
        assert_eq!(centromere_gap_lines(true).len(), 31);
    }
}
