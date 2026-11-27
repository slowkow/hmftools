//! ResultsWriter: the output files, filled gene group by gene group and written when the run ends.

use std::io::Write;

use crate::config::{LilacConfig, WriteType};
use crate::constants::*;

pub fn write_lines(path: &str, lines: &[String]) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    for l in lines {
        w.write_all(l.as_bytes())?;
        w.write_all(b"\n")?;
    }
    w.flush()
}

const REF_FRAGMENT_HEADER: &str = "Complex\tFragmentId\tFragmentCoords\tType\tFullAlleles\tWildAlleles";
const COMPLEX_HEADER: &str = "Genes\tScore\tComplexityPenalty\tComplexity\tHomozygousCount\tCohortFrequencyPenalty\tCohortFrequency\tRecoveryCount\tWildcardCount\tTotalCoverage\tUniqueCoverage\tSharedCoverage\tWildCoverage\tAlleleInfo";

/// One output file: its path and lines
pub struct OutFile {
    pub path: String,
    pub lines: Vec<String>,
}

impl OutFile {
    fn new(path: String, header: &str) -> std::io::Result<OutFile> {
        // created when the writer opens, as Java does; the content is written at the end
        std::fs::File::create(&path)?;
        Ok(OutFile { path, lines: vec![header.to_string()] })
    }
    fn write(&self) -> std::io::Result<()> {
        write_lines(&self.path, &self.lines)
    }
}

/// What one gene group adds to the output files, and the files it writes itself (in its order)
#[derive(Default)]
pub struct GroupOut {
    pub summary: Vec<String>,
    pub qc: Vec<String>,
    pub complexes: Vec<String>,
    pub fragments: Vec<String>,
    pub ref_fragments: Vec<String>,
    /// files written whole (path, content); a later gene group's file of the same name replaces it, as in Java
    pub files: Vec<(String, String)>,
}

pub struct Writers {
    pub summary: OutFile,
    pub qc: OutFile,
    pub complexes: OutFile,
    pub fragments: Option<OutFile>,
    pub ref_fragments: Option<OutFile>,
}

impl Writers {
    pub fn new(config: &LilacConfig) -> std::io::Result<Writers> {
        let frag = config.writes(WriteType::Fragments);
        // LilacAllele.generateFilename: basePath + File.separator + sample
        let summary = OutFile::new(format!("{}/{}.lilac.tsv", config.output_dir, config.sample), crate::qc::LILAC_ALLELE_HEADER)?;
        let qc = OutFile::new(format!("{}{}.lilac.qc.tsv", config.output_dir, config.sample), &crate::qc::LilacQC::header().join("\t"))?;
        // LI2: Java opens the candidate coverage file with the fragment file's header
        let complex_header = if crate::compat::java("LI2") { REF_FRAGMENT_HEADER } else { COMPLEX_HEADER };
        let complexes = OutFile::new(config.form_file_id(LILAC_FILE_CANDIDATE_COVERAGE), complex_header)?;
        let fragments = if frag {
            Some(OutFile::new(config.form_file_id(LILAC_FILE_FRAGMENTS), "Source\tReadId\tReadInfo\tGenes\tNucLociStart\tNucLociEnd\tAcidLociStart\tAcidLociEnd\tScope")?)
        } else {
            None
        };
        let ref_fragments = if frag { Some(OutFile::new(config.form_file_id(LILAC_FILE_CANDIDATE_FRAGS), REF_FRAGMENT_HEADER)?) } else { None };
        Ok(Writers { summary, qc, complexes, fragments, ref_fragments })
    }

    /// the lines of one gene group, appended in gene group order; its own files are written now
    pub fn append(&mut self, g: GroupOut) -> std::io::Result<()> {
        self.summary.lines.extend(g.summary);
        self.qc.lines.extend(g.qc);
        self.complexes.lines.extend(g.complexes);
        if let Some(f) = self.fragments.as_mut() {
            f.lines.extend(g.fragments);
        }
        if let Some(f) = self.ref_fragments.as_mut() {
            f.lines.extend(g.ref_fragments);
        }
        for (path, content) in g.files {
            std::fs::write(&path, content)?;
        }
        Ok(())
    }

    pub fn close(&self) -> std::io::Result<()> {
        if let Some(f) = &self.fragments {
            f.write()?;
        }
        if let Some(f) = &self.ref_fragments {
            f.write()?;
        }
        self.summary.write()?;
        self.qc.write()?;
        self.complexes.write()
    }
}

/// ResultsWriter.writeFragments: one line per fragment
pub fn fragment_line(source: &str, f: &crate::fragment::Fragment) -> String {
    let genes: Vec<String> = f.genes.iter().map(|g| g.long_name()).collect();
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        source,
        f.id(),
        f.read_info(),
        genes.join(";"),
        f.min_nucleotide_locus(),
        f.max_nucleotide_locus(),
        f.min_aa_locus(),
        f.max_aa_locus(),
        f.scope.name()
    )
}
