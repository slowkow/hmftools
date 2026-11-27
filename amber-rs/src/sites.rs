//! The AMBER sites (AmberSitesFile, as a TSV or a VCF), the -specific_regions and -target_regions_bed
//! filters of AmberApplication.loadAmberSites, and the BED resources.

use anyhow::{bail, Context, Result};
use std::io::BufRead;

use crate::config::{include_position, java_parse_double, java_parse_int, AmberConfig};
use crate::constants::TARGET_REGION_SITE_BUFFER;
use crate::types::{AmberBase, AmberSite};

/// A chromosome name, interned: the few distinct names of a run are kept for its lifetime.
pub fn intern(name: &str) -> &'static str {
    use std::sync::Mutex;
    static NAMES: Mutex<Vec<&'static str>> = Mutex::new(Vec::new());
    let mut names = NAMES.lock().unwrap();
    if let Some(n) = names.iter().find(|n| **n == name) {
        return n;
    }
    let n: &'static str = Box::leak(name.to_string().into_boxed_str());
    names.push(n);
    n
}

/// RefGenomeFunctions.stripChrPrefix
pub fn strip_chr(chromosome: &str) -> &str {
    chromosome.strip_prefix("chr").unwrap_or(chromosome)
}

/// HumanChromosome.contains: 1-22, X or Y, with or without the chr prefix
pub fn is_human_chromosome(chromosome: &str) -> bool {
    let t = strip_chr(chromosome);
    if t.bytes().all(|b| b.is_ascii_digit()) {
        // Java parses an empty string here and throws; treated as not human
        return matches!(t.parse::<i32>(), Ok(n) if (1..=22).contains(&n));
    }
    t == "X" || t == "Y"
}

/// The HumanChromosome of a human chromosome name (1-24 for 1-22, X, Y), as HumanChromosome.fromString;
/// 0 for anything else. The `chr` prefix is matched without regard to case, as fromString does.
pub fn human_chromosome(chromosome: &str) -> usize {
    let t = if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome };
    match t {
        "X" => 23,
        "Y" => 24,
        _ => match t.parse::<usize>() {
            Ok(n) if (1..=22).contains(&n) && !t.starts_with('0') && !t.starts_with('+') => n,
            _ => 0,
        },
    }
}

/// HumanChromosome.toString: the name without a prefix
pub fn human_chromosome_name(index: usize) -> &'static str {
    const NAMES: [&str; 25] = [
        "", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",
        "Y",
    ];
    NAMES[index]
}

/// RefGenomeVersion.versionedChromosome of a HumanChromosome
pub fn versioned_chromosome(v38: bool, index: usize) -> String {
    if v38 { format!("chr{}", human_chromosome_name(index)) } else { human_chromosome_name(index).to_string() }
}

/// Values keyed by HumanChromosome, as the Java ListMultimap<Chromosome, T>: one list per chromosome.
/// The key order of the Java maps is that of HashMap over enum identity hashes, which does not reach the
/// output; here the keys go in chromosome order.
#[derive(Clone, Debug)]
pub struct ByChr<T> {
    pub lists: Vec<Vec<T>>,
}

impl<T> Default for ByChr<T> {
    fn default() -> Self {
        ByChr { lists: (0..25).map(|_| Vec::new()).collect() }
    }
}

impl<T> ByChr<T> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn put(&mut self, chr_index: usize, value: T) {
        self.lists[chr_index].push(value);
    }
    pub fn get(&self, chr_index: usize) -> &[T] {
        &self.lists[chr_index]
    }
    pub fn len(&self) -> usize {
        self.lists.iter().map(|l| l.len()).sum()
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// (chromosome index, values) of the chromosomes with values
    pub fn chromosomes(&self) -> impl Iterator<Item = (usize, &Vec<T>)> {
        self.lists.iter().enumerate().filter(|(_, l)| !l.is_empty())
    }
    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.lists.iter().flat_map(|l| l.iter())
    }
    pub fn map<U>(&self, mut f: impl FnMut(&T) -> U) -> ByChr<U> {
        ByChr { lists: self.lists.iter().map(|l| l.iter().map(&mut f).collect()).collect() }
    }
    /// Multimaps.filterEntries
    pub fn filter(&self, mut f: impl FnMut(&T) -> bool) -> ByChr<T>
    where
        T: Clone,
    {
        ByChr { lists: self.lists.iter().map(|l| l.iter().filter(|x| f(x)).cloned().collect()).collect() }
    }
}

/// One AMBER site: AmberSite with its bases checked (Java's PositionEvidence fails on any other base).
#[derive(Clone, Copy, Debug)]
pub struct Locus {
    pub chr: &'static str,
    pub pos: i32,
    pub ref_: AmberBase,
    pub alt: AmberBase,
    pub snp_check: bool,
    pub gnomad_frequency: f64,
}

impl Locus {
    pub fn to_amber_site(&self) -> AmberSite {
        AmberSite::with_frequency(self.chr, self.pos, self.ref_.as_str(), self.alt.as_str(), self.snp_check, self.gnomad_frequency)
    }
}

/// The loci file and the order its sites came in.
pub struct Sites {
    pub by_chr: ByChr<Locus>,
}

fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    if path.ends_with(".gz") {
        Ok(Box::new(std::io::BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(std::io::BufReader::with_capacity(1 << 20, file)))
    }
}

fn parse_base(s: &str, line: &str) -> Result<AmberBase> {
    match AmberBase::from_str(s) {
        Some(b) => Ok(b),
        None => bail!("invalid base {s} in site: {line}"),
    }
}

fn locus(chr: &str, pos: i32, ref_: &str, alt: &str, snp_check: bool, freq: f64, line: &str) -> Result<Locus> {
    let r = parse_base(ref_, line)?;
    let a = parse_base(alt, line)?;
    if r == a {
        bail!("site with the same ref and alt base: {line}");
    }
    Ok(Locus { chr: intern(chr), pos, ref_: r, alt: a, snp_check, gnomad_frequency: freq })
}

/// AmberSitesFile.sites: a VCF when the first line says so, otherwise a TSV with the columns Chromosome,
/// Position, Ref, Alt, SnpCheck and optionally Frequency (0 without the column; 0.5 from a VCF).
pub fn load_sites(path: &str) -> Result<Sites> {
    let mut reader = open_text(path)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let header = header.trim_end_matches(['\n', '\r']).to_string();
    let mut by_chr = ByChr::new();
    if header.contains("fileformat=VCF") {
        load_vcf(reader, &mut by_chr)?;
        return Ok(Sites { by_chr });
    }
    let fields: Vec<&str> = header.split('\t').collect();
    let index = |name: &str| fields.iter().rposition(|f| *f == name);
    let (Some(chr_i), Some(pos_i), Some(ref_i), Some(alt_i), Some(snp_i)) =
        (index("Chromosome"), index("Position"), index("Ref"), index("Alt"), index("SnpCheck"))
    else {
        bail!("loci file {path} lacks a required column");
    };
    let freq_i = index("Frequency");
    let mut line = String::new();
    let mut count = 0usize;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let l = line.strip_suffix('\n').unwrap_or(&line);
        let l = l.strip_suffix('\r').unwrap_or(l);
        let values: Vec<&str> = l.split('\t').collect();
        let get = |i: usize| values.get(i).copied().with_context(|| format!("short line in {path}: {l}"));
        let chr = get(chr_i)?;
        if !is_human_chromosome(chr) {
            continue;
        }
        let ci = human_chromosome(chr);
        let pos = java_parse_int(get(pos_i)?).with_context(|| format!("invalid position in {path}: {l}"))?;
        let snp_check = get(snp_i)?.eq_ignore_ascii_case("true");
        let freq = match freq_i {
            Some(i) => java_parse_double(get(i)?).with_context(|| format!("invalid frequency in {path}: {l}"))?,
            None => 0.0,
        };
        by_chr.put(ci, locus(chr, pos, get(ref_i)?, get(alt_i)?, snp_check, freq, l)?);
        count += 1;
    }
    crate::info!("loaded {count} Amber germline sites from {path}");
    Ok(Sites { by_chr })
}

/// AmberSitesFile.loadVcf: unfiltered records on human chromosomes, the first alt allele, SNPCHECK from INFO
fn load_vcf(reader: Box<dyn BufRead>, by_chr: &mut ByChr<Locus>) -> Result<()> {
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            bail!("invalid VCF line: {line}");
        }
        // VariantContext.isFiltered: a FILTER other than PASS and '.'
        if f[6] != "." && f[6] != "PASS" {
            continue;
        }
        if !is_human_chromosome(f[0]) {
            continue;
        }
        let pos = java_parse_int(f[1]).with_context(|| format!("invalid VCF position: {line}"))?;
        let alt = f[4].split(',').next().unwrap_or(".");
        let snp_check = f[7].split(';').any(|kv| kv.split('=').next() == Some("SNPCHECK"));
        by_chr.put(human_chromosome(f[0]), locus(f[0], pos, f[3], alt, snp_check, 0.5, &line)?);
    }
    Ok(())
}

/// BedFileReader.loadBedFileChrMap: regions by HumanChromosome, lines with "Chromosome" skipped, start + 1.
/// AM5: Java starts a new list each time the chromosome changes and puts it in the map over any earlier
/// one, so a chromosome whose lines are not together keeps only its last group. Corrected: the groups
/// are appended.
pub fn load_bed_chr_map(path: &str) -> Result<ByChr<(i32, i32)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut map: ByChr<(i32, i32)> = ByChr::new();
    let mut current: usize = usize::MAX;
    for line in text.lines() {
        if line.contains("Chromosome") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            bail!("invalid slice BED entry: {line}");
        }
        if !is_human_chromosome(f[0]) {
            continue;
        }
        let ci = human_chromosome(f[0]);
        let (Some(s), Some(e)) = (java_parse_int(f[1]), java_parse_int(f[2])) else { bail!("invalid BED entry: {line}") };
        if ci != current {
            current = ci;
            if crate::compat::java("AM5") {
                map.lists[ci].clear();
            }
        }
        map.put(ci, (s + 1, e));
    }
    Ok(map)
}

/// AmberUtils.loadBedFromResource: chromosome, start + 1, end
pub fn bed_resource(text: &str) -> Vec<(String, i32, i32)> {
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.splitn(3, '\t').collect();
            (f[0].to_string(), f[1].parse::<i32>().unwrap() + 1, f[2].trim_end().parse::<i32>().unwrap())
        })
        .collect()
}

pub const TUMOR_ONLY_EXCLUDED_SNP_38: &str = include_str!("../resources/tumorOnlyExcludedSnp.38.bed");
pub const ROH_EXCLUDED_37: &str = include_str!("../resources/rohExcluded.37.bed");
pub const ROH_EXCLUDED_38: &str = include_str!("../resources/rohExcluded.38.bed");
pub const CENTROMERES_37: &str = include_str!("../resources/centromeres.37.tsv");
pub const CENTROMERES_38: &str = include_str!("../resources/centromeres.38.tsv");

/// AmberApplication.loadAmberSites: the loci file, restricted to -specific_regions / -specific_chr and to
/// the -target_regions_bed regions widened by 300 bases.
pub fn load_amber_sites(config: &AmberConfig) -> Result<ByChr<Locus>> {
    let mut sites = load_sites(&config.baf_loci_path)?.by_chr;
    if config.specific_regions.has_filters() {
        sites = sites.filter(|s| include_position(&config.specific_regions, s.chr, s.pos));
    }
    let Some(bed) = &config.target_regions_bed else { return Ok(sites) };
    let regions = match load_bed_chr_map(bed) {
        Ok(r) => r,
        Err(e) => {
            crate::error!("failed to load target regions file(): {{{bed}}}: {e}");
            std::process::exit(1);
        }
    };
    let mut result: ByChr<Locus> = ByChr::new();
    for (ci, regions) in regions.chromosomes() {
        let mut region_index = 0usize;
        let mut current = regions[0];
        for site in sites.get(ci) {
            if site.pos < current.0 - TARGET_REGION_SITE_BUFFER {
                continue;
            }
            while site.pos > current.1 + TARGET_REGION_SITE_BUFFER {
                region_index += 1;
                if region_index >= regions.len() {
                    break;
                }
                current = regions[region_index];
            }
            if region_index >= regions.len() {
                break;
            }
            if site.pos >= current.0 - TARGET_REGION_SITE_BUFFER && site.pos <= current.1 + TARGET_REGION_SITE_BUFFER {
                result.put(ci, *site);
            }
        }
    }
    Ok(result)
}

/// The sites of AmberApplication.hetLociTumorOnly: those outside the tumor-only excluded regions (38 only)
pub fn tumor_only_sites(config: &AmberConfig, sites: &ByChr<Locus>) -> ByChr<Locus> {
    let excluded = if config.v38 { bed_resource(TUMOR_ONLY_EXCLUDED_SNP_38) } else { Vec::new() };
    let mut removed = 0usize;
    let result = sites.filter(|s| {
        let out = excluded.iter().any(|(c, start, end)| c == s.chr && s.pos >= *start && s.pos <= *end);
        if out {
            removed += 1;
        }
        !out
    });
    if removed > 0 {
        crate::info!("removed {removed} blacklisted loci, {} remaining", result.len());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chromosomes() {
        assert!(is_human_chromosome("chr1") && is_human_chromosome("22") && is_human_chromosome("chrX"));
        assert!(!is_human_chromosome("chrM") && !is_human_chromosome("23") && !is_human_chromosome("chr1_alt") && !is_human_chromosome("chr"));
        assert_eq!(human_chromosome("chr10"), 10);
        assert_eq!(human_chromosome("Y"), 24);
        assert_eq!(versioned_chromosome(true, 23), "chrX");
        assert_eq!(versioned_chromosome(false, 2), "2");
        assert!(std::ptr::eq(intern("chr1"), intern("chr1")));
    }

    fn write(dir: &std::path::Path, name: &str, text: &str) -> String {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn loads_tsv_and_vcf() {
        let dir = std::env::temp_dir().join(format!("amber_rs_sites_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let tsv = write(
            &dir,
            "s.tsv",
            "Chromosome\tPosition\tRef\tAlt\tSnpCheck\tFrequency\nchr1\t100\tA\tG\tfalse\t0.25\nchrM\t5\tA\tG\tfalse\t0.1\nchr2\t50\tC\tT\tTRUE\t0.5\n",
        );
        let s = load_sites(&tsv).unwrap().by_chr;
        assert_eq!(s.len(), 2);
        assert_eq!(s.get(1)[0].gnomad_frequency, 0.25);
        assert!(s.get(2)[0].snp_check);
        let tsv2 = write(&dir, "s2.tsv", "Chromosome\tPosition\tRef\tAlt\tSnpCheck\n1\t100\tA\tG\tfalse\n");
        assert_eq!(load_sites(&tsv2).unwrap().by_chr.get(1)[0].gnomad_frequency, 0.0);
        let vcf = write(
            &dir,
            "s.vcf",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC,T\t.\tPASS\tSNPCHECK\nchr1\t20\t.\tA\tC\t.\tLOW\t.\nchr3\t30\t.\tG\tA\t.\t.\tX=1\n",
        );
        let v = load_sites(&vcf).unwrap().by_chr;
        assert_eq!(v.len(), 2);
        assert!(v.get(1)[0].snp_check && v.get(1)[0].alt == AmberBase::C && v.get(1)[0].gnomad_frequency == 0.5);
        assert!(!v.get(3)[0].snp_check);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn bed_groups() {
        let dir = std::env::temp_dir().join(format!("amber_rs_bed_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bed = write(&dir, "t.bed", "Chromosome\tStart\tEnd\n1\t10\t20\n2\t5\t6\n1\t30\t40\nchrUn\t1\t2\n");
        let fixed = load_bed_chr_map(&bed).unwrap();
        assert_eq!(fixed.get(1), &[(11, 20), (31, 40)]);
        let java = crate::compat::with_java(true, || load_bed_chr_map(&bed).unwrap());
        assert_eq!(java.get(1), &[(31, 40)]);
        assert_eq!(java.get(2), &[(6, 6)]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn resources() {
        let r = bed_resource(ROH_EXCLUDED_38);
        assert_eq!(r.len(), 8);
        assert_eq!(r[2], ("chr13".to_string(), 2, 17025624));
        assert_eq!(bed_resource(TUMOR_ONLY_EXCLUDED_SNP_38).len(), 33);
    }
}
