//! NeoEpitopeFinder (the jar's main class): NeoConfig, NeoSampleTask and the neo-epitope writer.
//!
//! For each sample: LINX's neo-epitope fusion candidates and PURPLE's somatic point mutations (missense,
//! frameshift and stop-lost impacts on a gene of the Ensembl cache) give candidate neo-epitopes, one per
//! qualifying transcript; per variant, candidates with the same or contained amino acids are merged and
//! the rest written to `<sample>.neo.neo_data.tsv`.

use anyhow::{bail, Context, Result};
use std::io::Write;

use hmf_common::config_builder::ConfigBuilder;
use hmf_common::ensembl::{EnsemblDataCache, TranscriptData};
use hmf_common::java_hash::JavaHashSet;
use hmf_common::neo::{
    check_add_dir_separator, neo_epitope_filename, neo_epitope_fusion_filename, neo_epitope_header, neo_epitope_line,
    read_neo_epitope_fusions, NeoEpitopeFusion, NeoEpitopeType,
};
use hmf_common::regions::{position_within, RefGenomeVersion};
use hmf_common::variant_context::VcfReader;
use hmf_common::{debug, error, info, warn};

use crate::compat;
use crate::epitope::*;

pub const DEFAULT_AMINO_ACID_REF_COUNT: i32 = 15;

/// NeoConfig.addCmdLineArgs
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_str("sample", false, "Sample - Id(s) separated by ';' or CSV file");
    cb.add_path("sample_id_file", false, "Sample ID CSV file with 'SampleId' column");
    cb.add_str("cancer_type", false, "Tumor cancer type (optional) - to retrieve cancer median TPM");
    cb.add_path("linx_dir", true, "Linx neoepitope directory");
    cb.add_path("somatic_vcf", true, "Purple somatic VCF (use '*') as required");
    cb.add_flag("write_trans_data", "Write transcript data for each neo-epitope");
    cb.add_integer("req_amino_acids", "Number of amino acids in neo-epitopes", DEFAULT_AMINO_ACID_REF_COUNT);
    cb.add_path("gene_id_file", false, "Restricted set of Gene IDs in CSV file");
    cb.add_str("ref_genome_version", true, "Ref genome version, 37 or 38");
    cb.add_path("ref_genome", true, "Path to reference genome fasta file");
    cb.add_path("ensembl_data_dir", false, "Ensembl data file directory");
    cb.add_logging_options();
    cb.add_str("output_dir", false, "Output directory");
    cb.add_str("output_id", false, "Output file ID");
    cb.add_thread_options(1);
}

/// NeoConfig
pub struct NeoConfig {
    pub samples: Vec<String>,
    pub ref_genome: FastaRef,
    pub ref_genome_version: RefGenomeVersion,
    pub restricted_gene_ids: Vec<String>,
    pub required_amino_acids: i32,
    pub write_trans_data: bool,
    pub linx_dir: Option<String>,
    pub somatic_vcf: Option<String>,
    pub output_dir: Option<String>,
    /// Java splits several samples between threads; each sample has its own file, so the port runs them in turn
    #[allow(dead_code)]
    pub threads: i32,
}

/// ConfigUtils.loadDelimitedIdFile: the id column of a CSV file, or its only column
pub fn load_delimited_id_file(filename: &str, id_column: &str) -> Vec<String> {
    let lines = match hmf_common::java_text::read_lines(filename) {
        Ok(l) => l,
        Err(e) => {
            error!("failed to read {} file({}): {}", id_column, filename, e);
            return Vec::new();
        }
    };
    let mut ids = Vec::new();
    let Some(first) = lines.first() else { return ids };
    let column_count = hmf_common::java_text::split_all(first, ',').len();
    let (index, skip) = if first.contains(id_column) {
        let f = hmf_common::java_text::fields_index_map(first, ',');
        match f.get(id_column) {
            Some(&i) => (i, 1),
            None => return ids,
        }
    } else {
        if column_count > 1 {
            return ids;
        }
        (0, 0)
    };
    for line in &lines[skip..] {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(v) = hmf_common::java_text::split_all(line, ',').get(index) {
            ids.push(v.to_string());
        }
    }
    ids
}

/// ConfigUtils.convertWildcardSamplePath
pub fn convert_wildcard_sample_path(path: &str, sample: &str) -> String {
    path.replace('*', sample)
}

impl NeoConfig {
    pub fn new(cb: &ConfigBuilder) -> Result<NeoConfig> {
        let mut samples = Vec::new();
        if let Some(s) = cb.value("sample") {
            samples.push(s);
        } else if let Some(f) = cb.value("sample_id_file") {
            samples = load_delimited_id_file(&f, "SampleId");
        }
        let ref_path = cb.value("ref_genome").context("missing ref_genome")?;
        let genome = hmf_common::refgenome::RefGenome::open(&ref_path).with_context(|| format!("failed to load ref genome({ref_path})"))?;
        let ref_genome_version = RefGenomeVersion::from_java(&cb.value("ref_genome_version").unwrap_or_default())?;
        let mut restricted_gene_ids = Vec::new();
        if let Some(f) = cb.value("gene_id_file") {
            restricted_gene_ids = load_delimited_id_file(&f, "GeneId");
        }
        Ok(NeoConfig {
            samples,
            ref_genome: FastaRef { genome },
            ref_genome_version,
            restricted_gene_ids,
            required_amino_acids: cb.integer("req_amino_acids"),
            write_trans_data: cb.has_flag("write_trans_data"),
            linx_dir: cb.value("linx_dir"),
            somatic_vcf: cb.value("somatic_vcf"),
            output_dir: cb.value("output_dir").map(|d| check_add_dir_separator(&d)),
            threads: cb.integer("threads"),
        })
    }
}

/// NeoEpitopeFinder.main after the command line is parsed
pub fn run(cb: &ConfigBuilder) -> Result<()> {
    let config = NeoConfig::new(cb)?;
    // NeoEpitopeFinder: the Ensembl cache with exons, every transcript, restricted gene ids
    let mut cache = EnsemblDataCache::new(&cb.value("ensembl_data_dir").unwrap_or_default(), config.ref_genome_version);
    cache.set_required_data(true, false, false, false);
    cache.set_restricted_gene_id_list(&config.restricted_gene_ids);
    if let Err(e) = cache.load(false) {
        // Java ignores the result of load(): the cache stays empty
        error!("{e:#}");
    }
    cache.create_gene_name_id_map();

    if config.samples.is_empty() {
        return Ok(());
    }
    if config.samples.len() == 1 {
        info!("processing sample({})", config.samples[0]);
    } else {
        info!("processing {} samples", config.samples.len());
    }
    for sample in &config.samples {
        let mut task = NeoSampleTask::new(sample, &config, &cache)?;
        task.process_sample()?;
    }
    info!("Neo-epitope annotations complete");
    Ok(())
}

/// A file whose text is kept in memory and written by `finish`: Java's BufferedWriter is closed only when
/// the sample completes, so a sample that fails leaves the file empty (for outputs below its 8 KB buffer)
pub struct DeferredFile {
    file: std::fs::File,
    buf: Vec<u8>,
}

impl DeferredFile {
    pub fn new(file: std::fs::File) -> DeferredFile {
        DeferredFile { file, buf: Vec::with_capacity(1 << 16) }
    }
    pub fn finish(mut self) -> std::io::Result<()> {
        self.file.write_all(&self.buf)
    }
}

impl Write for DeferredFile {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.buf.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// NeoSampleTask
pub struct NeoSampleTask<'a> {
    sample_id: String,
    config: &'a NeoConfig,
    cache: &'a EnsemblDataCache,
    next_neo_epitope_id: i32,
    writer: Option<DeferredFile>,
}

impl<'a> NeoSampleTask<'a> {
    pub fn new(sample_id: &str, config: &'a NeoConfig, cache: &'a EnsemblDataCache) -> Result<NeoSampleTask<'a>> {
        // NeoEpitopeFinder.initialiseNeoepitopeWriter
        let Some(output_dir) = &config.output_dir else {
            bail!("java.lang.NullPointerException: no output directory");
        };
        let writer = if output_dir.is_empty() {
            None
        } else {
            let filename = neo_epitope_filename(output_dir, sample_id);
            match std::fs::File::create(&filename) {
                Ok(f) => {
                    let mut w = DeferredFile::new(f);
                    writeln!(w, "{}", neo_epitope_header())?;
                    Some(w)
                }
                Err(e) => {
                    error!("error initialising neo-epitope output file: java.io.FileNotFoundException: {filename} ({e})");
                    None
                }
            }
        };
        Ok(NeoSampleTask { sample_id: sample_id.to_string(), config, cache, next_neo_epitope_id: 0, writer })
    }

    pub fn process_sample(&mut self) -> Result<()> {
        let fusions = self.get_sv_fusions();
        let point_mutations = self.get_somatic_variants()?;
        info!("sample({}) loaded {} fusions and {} point mutations", self.sample_id, fusions.len(), point_mutations.len());
        self.add_sv_fusions(&fusions)?;
        self.add_point_mutations(&point_mutations)?;
        if let Some(w) = self.writer.take() {
            w.finish()?;
        }
        Ok(())
    }

    fn get_somatic_variants(&self) -> Result<Vec<PointMutationData>> {
        let mut out = Vec::new();
        let Some(vcf) = &self.config.somatic_vcf else { return Ok(out) };
        let vcf = convert_wildcard_sample_path(vcf, &self.sample_id);
        if !std::path::Path::new(&vcf).exists() {
            warn!("Purple somatic VCF file({}) not found", vcf);
            return Ok(out);
        }
        let mut reader = VcfReader::open(&vcf).with_context(|| format!("invalid VCF file({vcf})"))?;
        while let Some(rec) = reader.next_record().with_context(|| format!(" failed to read somatic VCF file({vcf})"))? {
            if rec.is_filtered() {
                continue;
            }
            // must have a gene impact
            let mut valid = false;
            let mut gene_name = String::new();
            // VariantImpactSerialiser.fromAttributeValues: 10 values or an empty impact
            let impact = rec.info.as_string_list("IMPACT", "");
            if impact.len() == 10 {
                check_coding_effect(&impact[3])?;
                check_coding_effect(&impact[8])?;
                if hmf_common::java_text::parse_int(&impact[9]).is_none() {
                    bail!("java.lang.NumberFormatException: For input string: \"{}\"", impact[9]);
                }
                let (impact_gene, canonical_effect, worst_coding_effect) = (&impact[0], &impact[2], &impact[8]);
                if !impact_gene.is_empty()
                    && self.cache.gene_data_by_name(impact_gene).is_some()
                    && is_relevant_mutation(worst_coding_effect, canonical_effect)
                {
                    gene_name = impact_gene.clone();
                    valid = true;
                }
            }
            if !valid && rec.info.has("PAVE_TI") {
                // VariantTranscriptImpact.fromVariantContext: every impact is parsed first
                let mut impacts: Vec<(String, String)> = Vec::new();
                for ti in rec.info.as_string_list("PAVE_TI", "") {
                    let items: Vec<&str> = ti.split('|').collect();
                    if items.len() < 7 {
                        bail!("java.lang.ArrayIndexOutOfBoundsException: transcript impact({ti})");
                    }
                    if items.len() >= 10 && (hmf_common::java_text::parse_int(items[8]).is_none() || hmf_common::java_text::parse_int(items[9]).is_none()) {
                        bail!("java.lang.NumberFormatException: transcript impact({ti})");
                    }
                    impacts.push((items[1].to_string(), items[3].to_string()));
                }
                for (impact_gene, effects) in impacts {
                    if self.cache.gene_data_by_name(&impact_gene).is_none() {
                        continue;
                    }
                    if check_variant_effects(&effects) != "NONE" {
                        gene_name = impact_gene;
                        valid = true;
                        break;
                    }
                }
            }
            if !valid {
                continue;
            }
            out.push(PointMutationData {
                chromosome: rec.chrom.clone(),
                position: rec.pos,
                ref_bases: rec.ref_allele.clone(),
                alt: rec.alts.join(","),
                gene: gene_name,
                variant_copy_number: rec.info.as_double("PURPLE_VCN", 0.0),
                copy_number: rec.info.as_double("PURPLE_CN", 0.0),
                subclonal_likelihood: rec.info.as_double("SUBCL", 0.0),
            });
        }
        debug!("loaded {} somatic variants from file({})", out.len(), vcf);
        Ok(out)
    }

    fn get_sv_fusions(&self) -> Vec<NeoEpitopeFusion> {
        let Some(linx_dir) = &self.config.linx_dir else { return Vec::new() };
        let linx_dir = convert_wildcard_sample_path(linx_dir, &self.sample_id);
        let filename = neo_epitope_fusion_filename(&linx_dir, &self.sample_id);
        if !std::path::Path::new(&filename).exists() {
            warn!("Linx neo-epitope file({}) not found", filename);
            return Vec::new();
        }
        match read_neo_epitope_fusions(&filename) {
            Ok(f) => {
                debug!("loaded {} Linx neo-epitope candidates from file: {}", f.len(), filename);
                f
            }
            Err(_) => {
                error!("failed to read Linx neo-epitope file({})", filename);
                Vec::new()
            }
        }
    }

    /// the transcripts of a gene; Java's getTranscripts gives null for a gene not in the cache (outside
    /// -gene_id_file, say), and its streams then throw a NullPointerException (NE3)
    fn gene_transcripts(&self, gene_id: &str) -> Result<&'a [TranscriptData]> {
        match self.cache.transcripts(gene_id) {
            Some(l) => Ok(l.as_slice()),
            None if compat::java("NE3") => bail!("java.lang.NullPointerException: no transcripts for gene({gene_id})"),
            None => Ok(&[]),
        }
    }

    fn add_sv_fusions(&mut self, fusions: &[NeoEpitopeFusion]) -> Result<()> {
        for fusion in fusions {
            let mut list: Vec<NeoEpitope> = Vec::new();
            // find all transcripts where the breakend is inside the coding region on the 5' gene
            let up_list = self.cache.transcripts(&fusion.gene_ids[FS_UP]);
            let down_list = self.cache.transcripts(&fusion.gene_ids[FS_DOWN]);
            let up_names: Vec<&str> = fusion.transcripts[FS_UP].split(';').collect();
            let down_names: Vec<&str> = fusion.transcripts[FS_DOWN].split(';').collect();
            let same_gene = fusion.gene_ids[FS_UP] == fusion.gene_ids[FS_DOWN];

            // for same gene fusions, check that the SV isn't non-disruptive within the coding region of any other transcript(s)
            if same_gene && fusion.orientations[FS_UP] == ORIENT_FWD && fusion.orientations[FS_DOWN] == ORIENT_REV {
                let up = self.gene_transcripts(&fusion.gene_ids[FS_UP])?;
                if up.iter().any(|x| sv_is_non_disruptive_in_coding_transcript(&fusion.positions, x)) {
                    continue;
                }
            }
            for up_name in &up_names {
                let up_all = match up_list {
                    Some(l) => l.as_slice(),
                    None => self.gene_transcripts(&fusion.gene_ids[FS_UP])?,
                };
                let Some(up_trans) = up_all.iter().find(|x| x.trans_name == *up_name) else { continue };
                // NE4: Java unboxes the coding start and end of a non-coding upstream transcript: NullPointerException
                let (Some(cs), Some(ce)) = (up_trans.coding_start, up_trans.coding_end) else {
                    if compat::java("NE4") {
                        bail!("java.lang.NullPointerException: non-coding upstream transcript({})", up_trans.trans_name);
                    }
                    continue;
                };
                if !position_within(fusion.positions[FS_UP], cs, ce) {
                    continue;
                }
                for down_name in &down_names {
                    let down_all = match down_list {
                        Some(l) => l.as_slice(),
                        None => self.gene_transcripts(&fusion.gene_ids[FS_DOWN])?,
                    };
                    let Some(down_trans) = down_all.iter().find(|x| x.trans_name == *down_name) else { continue };
                    // if the same gene then must be the same transcript
                    if same_gene && up_trans.trans_id != down_trans.trans_id {
                        continue;
                    }
                    // must have a splice acceptor
                    if down_trans.exons.len() <= 1 {
                        continue;
                    }
                    let (range_start, range_end) = if down_trans.strand == POS_STRAND {
                        (down_trans.trans_start - DOWNSTREAM_PRE_GENE_DISTANCE, down_trans.trans_end)
                    } else {
                        (down_trans.trans_start, down_trans.trans_end + DOWNSTREAM_PRE_GENE_DISTANCE)
                    };
                    if !position_within(fusion.positions[FS_DOWN], range_start, range_end) {
                        continue;
                    }
                    let mut ne = NeoEpitope::fusion(fusion);
                    ne.set_transcript_data(up_trans, down_trans);
                    list.push(ne);
                }
            }
            self.process_neo_epitopes(list)?;
        }
        Ok(())
    }

    fn add_point_mutations(&mut self, point_mutations: &[PointMutationData]) -> Result<()> {
        for pm in point_mutations {
            let gene = self.cache.gene_data_by_name(&pm.gene).context("gene data")?;
            let mut list: Vec<NeoEpitope> = Vec::new();
            for trans in self.gene_transcripts(&gene.gene_id)? {
                let (Some(cs), Some(ce)) = (trans.coding_start, trans.coding_end) else { continue };
                if !position_within(pm.position, cs, ce) {
                    continue;
                }
                // ignore IG and TR transcripts
                if is_immune_transcript(trans) {
                    continue;
                }
                // check for a mutation within the stop codon at the bounds of the transcript
                if pm.position <= trans.trans_start || pm.position >= trans.trans_end {
                    continue;
                }
                // must be exonic
                if !trans.exons.iter().any(|x| position_within(pm.position, x.start, x.end)) {
                    continue;
                }
                if trans.bio_type == hmf_common::ensembl::BIOTYPE_NONSENSE_MED_DECAY {
                    continue;
                }
                let mut ne = NeoEpitope::point_mutation(pm);
                ne.set_transcript_data(trans, trans);
                list.push(ne);
            }
            self.process_neo_epitopes(list)?;
        }
        Ok(())
    }

    fn process_neo_epitopes<'b>(&mut self, mut list: Vec<NeoEpitope<'b>>) -> Result<()> {
        if list.is_empty() {
            return Ok(());
        }
        let ref_genome = &self.config.ref_genome;
        let req = self.config.required_amino_acids;
        for ne in list.iter_mut() {
            ne.set_coding_bases(ref_genome, req)?;
        }
        for ne in list.iter_mut() {
            ne.set_amino_acids(ref_genome, req)?;
        }
        for ne in list.iter_mut() {
            ne.set_nonsense_mediated_decay();
        }
        for ne in list.iter_mut() {
            ne.set_skipped_splice_sites(self.cache);
        }

        // consolidate duplicates
        let mut i = 0;
        while i < list.len() {
            {
                let ne = &list[i];
                // filter out NEs with a novel stop codon
                if ne.novel_acid == STOP_SYMBOL {
                    i += 1;
                    continue;
                }
                if !ne.valid {
                    debug!("skipping invalid neo: {}", ne.describe());
                    i += 1;
                    continue;
                }
                // filter out missense if in some transcripts their NE is the same as the wild-type
                if ne.variant_type() == NeoEpitopeType::Missense
                    && !ne.wildtype_acids().is_empty()
                    && ne.amino_acid_string().contains(ne.wildtype_acids())
                {
                    i += 1;
                    continue;
                }
            }
            let mut up_names: JavaHashSet<String> = JavaHashSet::new();
            let mut down_names: JavaHashSet<String> = JavaHashSet::new();
            up_names.insert(list[i].trans(FS_UP).trans_name.clone());
            down_names.insert(list[i].trans(FS_DOWN).trans_name.clone());
            let aa = list[i].amino_acid_string();

            if !self.config.write_trans_data {
                let mut j = i + 1;
                while j < list.len() {
                    let (a, b) = (&list[i], &list[j]);
                    let min_nmd = a.nmd_bases_min.min(b.nmd_bases_min);
                    let max_nmd = a.nmd_bases_max.max(b.nmd_bases_max);
                    let min_cbl = a.coding_bases_length_min.min(b.coding_bases_length_min);
                    let max_cbl = a.coding_bases_length_max.max(b.coding_bases_length_max);
                    let other_aa = b.amino_acid_string();
                    // remove exact matches or take the longer if one is a subset
                    let other = if aa.contains(&other_aa) {
                        list.remove(j)
                    } else if other_aa.contains(&aa) {
                        let other = list.remove(j);
                        list[i] = other.clone();
                        other
                    } else {
                        j += 1;
                        continue;
                    };
                    // take the shortest NMD base count
                    let ne = &mut list[i];
                    ne.nmd_bases_min = min_nmd;
                    ne.nmd_bases_max = max_nmd;
                    ne.coding_bases_length_min = min_cbl;
                    ne.coding_bases_length_max = max_cbl;
                    // collect up all transcripts
                    up_names.insert(other.trans(FS_UP).trans_name.clone());
                    down_names.insert(other.trans(FS_DOWN).trans_name.clone());
                }
            }
            let ne_id = self.next_neo_epitope_id;
            self.next_neo_epitope_id += 1;
            if let Some(w) = self.writer.as_mut() {
                let file = list[i].to_file(ne_id, &up_names, &down_names)?;
                writeln!(w, "{}", neo_epitope_line(&file))?;
            }
            i += 1;
        }
        Ok(())
    }
}

/// CodingEffect.valueOf: IllegalArgumentException for an unknown name
fn check_coding_effect(s: &str) -> Result<()> {
    match s {
        "NONSENSE_OR_FRAMESHIFT" | "SPLICE" | "MISSENSE" | "SYNONYMOUS" | "NONE" | "UNDEFINED" => Ok(()),
        _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.CodingEffect.{s}"),
    }
}
