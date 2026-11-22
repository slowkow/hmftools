//! RefGenomeInterface: reference bases by 1-based inclusive coordinates.

use std::sync::Arc;

use crate::common::refgenome::RefGenome;

pub trait RefGenomeInterface: Send + Sync {
    /// null in Java (None) when the start is below 1 or the end is before the start
    fn get_bases(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> Option<Vec<u8>>;
    fn get_base_string(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> Option<String> {
        self.get_bases(chromosome, pos_start, pos_end).map(|b| b.iter().map(|c| *c as char).collect())
    }
    fn get_chromosome_length(&self, chromosome: &str) -> i32;
    /// contigs in dictionary order with their lengths
    fn chromosomes(&self) -> Vec<(String, i32)>;
}

pub type RefGenomeRef = Arc<dyn RefGenomeInterface>;

/// The indexed FASTA (Java RefGenomeSource).
pub struct RefGenomeSource {
    /// the reader is used by one thread at a time
    genome: std::sync::Mutex<RefGenome>,
    chromosomes: Vec<(String, i32)>,
}

unsafe impl Send for RefGenomeSource {}
unsafe impl Sync for RefGenomeSource {}

impl RefGenomeSource {
    pub fn open(path: &str) -> anyhow::Result<RefGenomeSource> {
        let genome = RefGenome::open(path)?;
        let chromosomes = genome.sequence_names().iter().map(|n| (n.clone(), genome.chromosome_length(n))).collect();
        Ok(RefGenomeSource { genome: std::sync::Mutex::new(genome), chromosomes })
    }
}

impl RefGenomeInterface for RefGenomeSource {
    fn get_bases(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> Option<Vec<u8>> {
        if pos_start < 1 || pos_end < pos_start {
            return None;
        }
        Some(self.genome.lock().unwrap().bases(chromosome, pos_start, pos_end).unwrap_or_else(|e| panic!("reference read failed for {chromosome}:{pos_start}-{pos_end}: {e}")))
    }
    fn get_chromosome_length(&self, chromosome: &str) -> i32 {
        self.chromosomes.iter().find(|(n, _)| n == chromosome).map(|(_, l)| *l).unwrap_or(0)
    }
    fn chromosomes(&self) -> Vec<(String, i32)> {
        self.chromosomes.clone()
    }
}
