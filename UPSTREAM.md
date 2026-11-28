# Following upstream hmftools

The Rust ports reproduce hmftools as it was at a known commit of
[hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools). Upstream keeps changing. We follow it
on purpose, one reviewed change at a time, rather than by merging whatever arrives: a change is taken
only when it is sensible and we can defend it, and a fault that arrives with a change is not copied into
the port's default behaviour.

This file records where each port started, the rules for following upstream, and every upstream
change since, with what we decided and why.

## Where the ports started

| Port | Java program | Started from upstream commit | Date of that commit |
|---|---|---|---|
| `sage-rs` | SAGE 5.1 | `a224a2cd8e7797673449c2ef1cbca19b6007fe01` (tag `sage-v5.1-beta.4`) | 2026-11-17 |
| `redux-rs` | REDUX 2.1 | `a224a2cd8e7797673449c2ef1cbca19b6007fe01` | 2026-11-17 |
| `esvee-rs` | ESVEE 2.1 | `a224a2cd8e`; Java reference rebuilt at `ca04e6c39d`, which adds one ESVEE fix (`76f6774e65`) | 2026-11-20 |
| every other port | AMBER 4.4, COBALT 3.1, BamTools 1.6.2, TEAL 1.4, LILAC 2.1, CIDER 1.2, PAVE 2.0, PURPLE 4.5, LINX 2.4, SIGS 1.3, Virus Interpreter 1.7.1, PEACH 2.0, CHORD 2.1.2, NEO 1.4, CUPPA 2.5 (with pycuppa), qsee 1.0, ORANGE 5.1, ISOFOX 2.1 (and 2.0.1) | `ca04e6c39d9a320ee3f65b5f297bbb10adb4cc71` | 2026-11-20 |
| `orange-rs` | ORANGE 5.1 | `ca04e6c39d`, plus upstream's Viridian loading (`f5f660cc9f`), taken while porting: ORANGE as at `883ab11bba` | 2026-11-24 |
| `hmf-common-rs` | the parts of hmf-common the ports use | `ca04e6c39d` | 2026-11-20 |
| `hmf-plot-rs` | draws what PURPLE's, LINX's, CUPPA's and qsee's R, Circos and matplotlib code draws, and ORANGE's PDF | the scripts at `ca04e6c39d`; BUI pinned at `1af5d5084d92e8ed30c637076bc964c3769386ca` | 2026-11-20 |
| `markdups-rs` | Picard MarkDuplicates, not an hmftools program | Picard 3.3.0 with htsjdk 4.1.3, as bundled in GATK 4.6.1.0 (the image oncoanalyser runs) | (Picard release) |

The `rust-ports` branch was cut from `a224a2cd8e` (first commit `a0b24b3299`, SAGE, 2026-11-18). Upstream
master at `ca04e6c39d` was merged into it once, as `4184ab9cb4` on 2026-11-22; the Java sources on the
branch have been that commit's since, and the Java jars every later port was checked against were built
from them.

ISOFOX's 2.0.1 mode follows the released ISOFOX 2.0.1 (the image `hmftools-isofox:2.0.1--hdfd78af_0`), which
oncoanalyser runs on STAR alignments.

## What each port tracks now

The commit whose Java `-java_compat` reproduces. A port not listed still tracks its start commit.

| Port | Tracks | Since | What moving brought |
|---|---|---|---|
| `sage-rs` | `883ab11bba` | 2026-11-28 | the scaled base-quality penalty, SAGE's own NM count, `##sageAppendVersion`, and the soft-clip realignment index with one correction (S57) |
| `redux-rs` | `883ab11bba` | 2026-11-28 | nothing (REDUX unchanged upstream; reviewed with SAGE) |
| `orange-rs` | `883ab11bba` | 2026-11-27 (while porting) | Viridian and Virus Interpreter loading |
| `purple-rs` | `883ab11bba` | 2026-11-27 | the new biallelic model, with one correction (PU10) |
| `linx-rs` | `883ab11bba` | 2026-11-28 | enhancer fusions up to 1 Mb, with one correction (LX6); exonic 5' breakends fuse more widely |
| `pave-rs` | `883ab11bba` | 2026-11-28 | the MET exon 14 special splice region |
| `qsee-rs` | `883ab11bba` | 2026-11-28 | exon coverage QC; stranded base-quality recalibration features and the two-page report, with one correction (QS17) |

## How we follow upstream

1. **What we track.** Upstream master. Hartwig also publishes pre-releases from side branches before they
   reach master; those are listed below when we see them and reviewed when they are on master, unless the
   pipeline needs them sooner.
2. **Every change is reviewed.** Each upstream commit that touches a ported program, or code in
   hmf-common that a port uses, is read in full (message, code, tests, documentation) and gets one
   decision per change it holds:
   - **accepted**: ported as upstream wrote it;
   - **accepted with a correction**: ported, with a fault in the new code corrected by default; the fault
     gets a finding id like any other Java fault, and `-java_compat` reproduces upstream;
   - **rejected**: not taken into the default, with the reason; `-java_compat` still reproduces upstream,
     so the port stays comparable with Java at the tracked commit;
   - **held**: the purpose or the effect is unclear; asked about before deciding;
   - **not applicable**: code no port runs (another program, a function oncoanalyser does not use,
     documentation only).
3. **What makes a change defensible.** Its purpose is stated or clear from the code and tests; it does
   what it says; it does not break the cases it was not aimed at; and its effect on our test cases is
   what the purpose predicts. For changes of a model or a threshold, the effect is measured on our cases
   (outputs before and after, in the default mode) and judged with whatever evidence there is: upstream's
   rationale and tests, the data, and a truth set where one exists. A rejection needs a concrete reason:
   a fault, a contradiction with the documentation, or worse results.
4. **The tracked commit.** A port's tracked commit moves only when every upstream change up to it has a
   decision. The Java reference for `-java_compat` is the jar built at the tracked commit, and the port's
   regression cases are remade with it. Old references are kept.
5. **Versions.** Upstream does not always change the version number when it changes behaviour (SAGE stayed
   5.1 from beta.4 to beta.5 and master), so the commit, not the version, says what a port reproduces.
6. **Pipeline comparisons.** The Java images oncoanalyser runs are pre-releases, some from side branches. A
   comparison inside the pipeline is only exact when the image and the port's tracked commit agree; where
   they do not, the comparison is made with a jar built at the tracked commit.

## Upstream changes since the start

Target of the current review: upstream master `883ab11bba24aa6d3bc28317dfe953331b583445` (2026-11-27). Java
jars at that commit, for the references, are built in a separate checkout.

### Since `a224a2cd8e` (SAGE, REDUX, ESVEE)

| Upstream commit | Date | Program | Change | Decision |
|---|---|---|---|---|
| `de50386378` (tag `sage-v5.1-beta.5`) | 2026-11-19 | SAGE | base-quality penalty scaled from 6 (recalibrated quality 24 and below) to 12 (36 and above), in place of a fixed 12 | **accepted**, on upstream's rationale: on our data it changes no call a truth set can judge. On the GIAB HG008-T somatic draft benchmark (v0.4), chr21 PASS calls are the same before and after (SBX 141 true, 8 false, 38 missed; Ultima 103, 19, 76); Illumina chr21 gains one PASS call (6,023 to 6,385 records) |
| `de50386378` | 2026-11-19 | SAGE | SAGE counts NM itself for reads without an NM tag | **accepted**: it fixes the crash of sage-rs S54 in the evidence stage. The candidate stage still passes no reference genome, so S54's correction stays there |
| `de50386378` | 2026-11-19 | SAGE | append: the input VCF must come from SAGE 4.0 or later; `##sageAppendVersion` written | **accepted**; now written in both modes, so WISP applies the threshold meant for scaled-penalty qualities. Noted, kept as upstream: a second append leaves the first one's version line (S59); versions compared as decimal numbers, so "5.10" reads as 5.1 (S60) |
| `de50386378` | 2026-11-19 | SAGE | realigned read index for inserts before a soft clip | **accepted with a correction** (sage-rs S57): upstream moves the core end back by the whole insert length whenever the read context's CIGAR ends in a soft clip; that is right only when the insert itself was soft-clipped. With the insert held as an `I` element the shift loses reads that need realigning (9 and 7 inserts in normal and tumor on Illumina chr21, germline mode); with part of a long insert held as an `I` it shifts one base too far and counts reads of a neighbouring allele (14 tumor reads at one public SBX site). The port shifts by the inserted bases the CIGAR does not hold as an insert. Somatic runs with soft filters are unchanged |
| `de50386378` | 2026-11-19 | SAGE | README | noted (S58): it still documents a fixed penalty of 12 and a `base_qual_fixed_penalty` option the code never had |
| `de50386378` | 2026-11-19 | WISP | purity changes | not applicable (not ported) |
| `76f6774e65` | 2026-11-19 | ESVEE, hmf-common | prep BAM loading with `-output_id`; specific-region position check without chromosomes | ESVEE: accepted (in the ESVEE reference since `ca04e6c39d`). The hmf-common check: accepted, no change needed: every hmftools caller checks for filters first, and hmf-common-rs, sage-rs and redux-rs never had the old behaviour |
| `c2564fce30` | 2026-11-19 | qsee, BamTools | limit of detection for exon variants | accepted: before `ca04e6c39d`, so in qsee-rs and bam-tools-rs from the start |
| `a6a26719c4` | 2026-11-18 | Viridian | viral read alignments | not applicable (not ported) |

REDUX itself has not changed upstream since `a224a2cd8e`.

### Since `ca04e6c39d` (all other ports)

| Upstream commit | Date | Program | Change | Decision |
|---|---|---|---|---|
| `7a4490caa6` | 2026-11-23 | LINX | enhancer (IG and other) fusions: the 3' gene's upstream distance for enhancer partners 250 kb to 1 Mb | **accepted with a correction** (linx-rs LX6): Java applies the 1 Mb to every gene near a breakend in an enhancer partner's window, not only to the enhancer partner, so FGFR1-3 fusions, which skip the distance check, can take partners 300-740 kb away that no other LINX rule allows. The port gives the 1 Mb to the enhancer partner genes only. Evidence: IGH-MYC at 400 kb and IGL-BCL6 at 700 kb are now reported, nothing within the old limits changes; on synthetic FGFR2 cases Java reports far partners in 4 of 7 and replaces a 5 kb partner by a 736 kb one in another. On the real cases no fusion, driver or visualiser file changes; `breakend.tsv` loses the unrelated-gene rows the old 250 kb windows also gave (17 on HG008) |
| `7a4490caa6` | 2026-11-23 | LINX | three fusion settings fixed, and with them `-skip_unphased_fusions` removed | accepted (oncoanalyser does not pass the option; the port now refuses it, as Java does) |
| `7a4490caa6` | 2026-11-23 | LINX | a 5' breakend inside a coding exon may fuse with any valid 3' transcript, not only one on the same SV | accepted: consistent with how LINX treats such breakends; on a synthetic BCR-ABL1 it finds the canonical in-frame exon 13 to exon 2 fusion where it chose an out-of-frame non-coding transcript before |
| `7a4490caa6` | 2026-11-23 | LINX | visualiser data: the enhancer side of a fusion drawn over its enhancer region (all enhancers) | accepted: TRA and BCL11B fusions get their row. The pipeline's known fusion data already has the 53 enhancer pairs |
| `110f42ee7d` | 2026-11-23 | COMPAR | compares Viridian reports | not applicable (not ported; not in oncoanalyser) |
| `f5f660cc9f` | 2026-11-24 | ORANGE, hmf-common | loads Viridian and Virus Interpreter files | accepted in `orange-rs` (hmftools `e471239079`) |
| `3807790a20` | 2026-11-24 | PURPLE | README: purity status and QC column descriptions | not applicable (documentation). The corrections are right, checked against the files PURPLE writes; purple-rs's README did not repeat the old descriptions |
| `4adc014589` | 2026-11-24 | PAVE | special splice regions (MET exon 14, intronic bases -30 to -10), with a new resource file | **accepted**. The entry is right: MET's canonical transcript (ENST00000397752, + strand), exon 14 is the in-frame exon lost in METex14 skipping, and -30 to -10 before its acceptor is intron 13's polypyrimidine tract and branch point. The rule stays on that transcript and exon; only the splice-region flag changes, not effects or HGVS. The resource is built into the binary. Two edges are questions for upstream, kept as upstream (pave-rs PA10): base -9 is in neither the standard region (-1 to -8) nor the new one, and insertions and deletions count when only their flanking base touches the region |
| `91286d6867` | 2026-11-24 | COBALT | panel normalisation builder: a region is also filtered when its 10th percentile ratio is below 0.1 | not applicable (the builder makes a panel resource; oncoanalyser does not run it and it is not ported) |
| `bc9b187516` | 2026-11-24 | PURPLE | biallelic model for somatic variants: the probability of LOH grows with purity, the no-wildtype test allows for allele-fraction uncertainty and is damped at high depth, and the variant copy number itself counts as evidence | **accepted with a correction** (purple-rs PU10: the threshold read count wrapped when converted to an int, giving the wrong probability in segments with copy number 0 or below). Evidence: every biallelic flag it changes on HG008 (260 of 260 tumor/normal) agrees with the reads alone; on the private pair 7 agree, 25 are inconclusive, none disagrees; PASS variants whose flag contradicts the reads fall from 478 to 361 (HG008) and 10 to 5 (private pair). Noted and kept as upstream (PU11): the copy-number step ignores read depth, so a 5-read variant can score 1.0. Upstream's README still describes the old model. purple-rs tracks `883ab11bba` since 2026-11-27 |
| `bc9b187516` | 2026-11-24 | ORANGE, PEACH, SAGE, TAUR | pom names | not applicable |
| `bdceb2302a` | 2026-11-26 | qsee | exon coverage QC: `gene_lod` becomes `exon_coverage_qc` (limit-of-detection VAF floored at 0.01), exon statuses PASS / WARN / EXCLUDED from coverage alone, a "Low coverage exons" summary row, missed-variant likelihood from exon coverage, `-lod_min_alt_reads`, `-exon_exclusions_file` and an exclusions builder | **accepted**. It also fixes upstream's `-targeted_mode` inversion (qsee-rs QS2), so the port's switch for it is gone. On our data the summary row is WARN for both tumors (43% and 73% of exons) at upstream's thresholds, which assume about 100x tumors; missed-variant values are unchanged |
| `883ab11bba` | 2026-11-27 | qsee | read strand in the base-quality recalibration features and plots; the report on two pages | **accepted with a correction** (qsee-rs QS17): the R script now writes `.page_1.png` and `.page_2.png` and no longer the `.qsee.vis.report.png` that ORANGE at the same commit embeds, so ORANGE would show "Plot not available". The port also writes page 1 under the old name by default. The pipeline's qsee cohort file predates the new feature names, so the report shows no cohort boxes for them until the cohort file is remade |

### Released from side branches, not on master yet (2026-11-27)

| Release tag | Commit | Branch | Changes beyond master | Decision |
|---|---|---|---|---|
| `purple-v4.5-beta.6` | `2ed3ac0401` | `AUS522-purple-minors` | TMB-high threshold 10 to 13 (settable as `tmb_high_threshold`); tumor-only TMB from all variants; hotspot SVs in tumor-only mode; one QC warning skipped in targeted mode | waits for master. Reviewed ahead: `c5396458be` the option yes, the default of 13 held (no rationale given; the branch README still says 10; no effect on our cases); `8f895807a6` (no copy-number-noise warning for panels) yes; `e2588c9fee` (tumor-only TMB with the panel rules) held: it drops indels from tumor-only TMB though tumor/normal TMB counts them, about -40% on our tumor-only case, and removes `-target_coding_factor`; its int-division fix is right on its own; `2ed3ac0401` (hotspot SVs in the somatic-only fit) yes |
| `orange-v5.1-beta.3` | `24642c2e1e` | `AUS482-orange-minors` | RNA sample id in the JSON; signatures in the JSON limited to the existing thresholds; unused SAGE germline option removed | waits for master |
| `esvee-v2.1-beta.4` | `261e6f9f98` | (deleted) | the same fix as `76f6774e65` | covered |
