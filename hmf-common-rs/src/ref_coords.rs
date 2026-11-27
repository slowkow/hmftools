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
}
