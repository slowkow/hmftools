import java.awt.Color;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.*;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.vis.*;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.vis.AminoAcidEvent;
import com.hartwig.hmftools.sage.vis.AminoAcidUtil;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;

// Prints the expected values of sage-rs/src/tests/visualiser.rs, one "key<TAB>value" per line.
// Not part of any build. To regenerate vis_java_fixtures.tsv.gz (the argument is where the synthetic
// FASTA is written; log lines, which start with a time, are dropped):
//   javac -cp sage/target/sage-5.1-jar-with-dependencies.jar -d /tmp/visfix sage-rs/src/tests/data/VisFixtures.java
//   java -cp sage/target/sage-5.1-jar-with-dependencies.jar:/tmp/visfix VisFixtures /tmp/visfix/synthetic.fasta \
//       | grep -v '^[0-9][0-9]:' | gzip -9 -n > sage-rs/src/tests/data/vis_java_fixtures.tsv.gz
// The font metrics, and the text positions in the SVG that follow from them, are those of the machine it runs on.
public class VisFixtures
{
    static void out(String key, Object value) { System.out.println(key + "\t" + value); }

    static String genome(int length)
    {
        StringBuilder sb = new StringBuilder();
        long state = 12345;
        for(int i = 0; i < length; i++)
        {
            state = (state * 1103515245L + 12345L) & 0x7fffffffL;
            sb.append("ACGT".charAt((int) ((state >> 16) & 3)));
        }
        return sb.toString();
    }

    static String models(List<GeneRegionViewModel> models)
    {
        StringJoiner sj = new StringJoiner(" ");
        for(GeneRegionViewModel m : models)
        {
            if(m instanceof GeneRegionViewModel.AminoAcidViewModel a)
                sj.add("A:" + a.region().start() + "-" + a.region().end() + ":" + a.aminoAcidPos() + ":" + a.ref() + ":" + a.alt() + ":" + a.basesBeforeRightInsert() + ":" + a.isStart());
            else if(m instanceof GeneRegionViewModel.NonCodingExonicRegionViewModel n)
                sj.add("N:" + n.region().start() + "-" + n.region().end() + ":" + n.basesBeforeRightInsert());
            else if(m instanceof GeneRegionViewModel.IntronicRegionViewModel n)
                sj.add("I:" + n.region().start() + "-" + n.region().end() + ":" + n.basesBeforeRightInsert());
            else
                sj.add("D:" + m.region().start() + "-" + m.region().end());
        }
        return sj.toString();
    }

    static String events(List<AminoAcidEvent> events)
    {
        StringJoiner sj = new StringJoiner(" ");
        for(AminoAcidEvent e : events)
        {
            if(e instanceof AminoAcidEvent.AminoAcidSubstitution s) sj.add("S:" + s.aminoAcidStartPos() + ":" + s.alt());
            else if(e instanceof AminoAcidEvent.AminoAcidIns s) sj.add("I:" + s.aminoAcidStartPos() + ":" + s.length());
            else if(e instanceof AminoAcidEvent.AminoAcidDel s) sj.add("D:" + s.aminoAcidStartPos() + ":" + s.aminoAcidEndPos());
            else if(e instanceof AminoAcidEvent.AminoAcidExt s) sj.add("E:" + s.aminoAcidStartPos());
            else sj.add("L:" + e.aminoAcidStartPos());
        }
        return sj.toString();
    }

    static TranscriptData transcript(byte strand, int codingStart, int codingEnd)
    {
        TranscriptData t = new TranscriptData(1, "ENST0001", "ENSG0001", true, strand, 1100, 1300, codingStart, codingEnd, "protein_coding", "");
        List<ExonData> exons = new ArrayList<>();
        exons.add(new ExonData(1, 1100, 1130, strand == 1 ? 1 : 3, -1, -1));
        exons.add(new ExonData(1, 1160, 1210, 2, -1, -1));
        exons.add(new ExonData(1, 1240, 1300, strand == 1 ? 3 : 1, -1, -1));
        t.setExons(exons);
        return t;
    }

    public static void main(String[] args) throws Exception
    {
        // CSS
        CssBuilder e = CssBuilder.EMPTY;
        out("css.empty", "[" + e + "]");
        CssBuilder border = e.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        out("css.table", border.marginRight(CssSize.px(10.0)));
        CssBuilder header = e.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        out("css.merge", e.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(border).merge(header));
        out("css.hidden", e.padding(CssSize.ZERO).margin(CssSize.ZERO).display("none"));
        out("css.shadow", e.fontSizePt(8).height(CssSize.em(1)).leftBoxShadow(CssSize.px(2), Color.RED).rightBoxShadow(CssSize.px(1.5), Color.BLUE));
        out("css.shadow_only", e.leftBoxShadow(CssSize.px(2), Color.RED));
        out("css.pct", e.fontSizePct(0.123456));
        out("css.override", e.textAlign("left").color(ColorUtil.PURPLE).merge(e.textAlign("right")));
        out("css.misc", e.borderBottom(CssSize.px(2), "solid", Color.BLACK).verticalAlign("top").writingMode("vertical-rl").borderSpacing(CssSize.ZERO).noBorder().fontFamily("sans-serif"));
        out("size.tiny", CssSize.px(0.00004));
        out("size.small", CssSize.px(0.00005));
        out("size.em", CssSize.em(1.23456));
        out("size.neg", CssSize.px(2).scale(-1));
        out("size.half", CssSize.em(0.5f));

        // font metrics
        for(int size : new int[] { 18, 26, 30 })
        {
            java.awt.Font font = new java.awt.Font("SansSerif", java.awt.Font.PLAIN, size);
            for(String s : new String[] { "7,920,552", "G", "12", "*", "W", "" })
            {
                java.awt.geom.Rectangle2D r = SvgUtil.getStringBounds(font, s);
                out("font." + size + "." + s, r.getWidth() + " " + r.getHeight() + " " + SvgUtil.getFontDescent(font));
            }
        }

        // colours
        Color c = ColorUtil.interpolateColors(Color.RED, Color.GREEN, 37 / 60.0);
        out("color.interpolate", c.getRed() + "," + c.getGreen() + "," + c.getBlue());
        c = ColorUtil.lighten(new Color(151, 83, 9), 1.0 - 1.0 * 11 / 37);
        out("color.lighten", c.getRed() + "," + c.getGreen() + "," + c.getBlue());
        c = SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR;
        out("color.overlap", c.getRed() + "," + c.getGreen() + "," + c.getBlue());

        // SVG
        out("svg.box", SvgRender.renderColoredBox(12.0, new Color(1, 2, 3)).getSVGElement());
        out("svg.coords", SvgRender.renderCoords(12.0, new BaseRegion(995, 1012), 1003, 10).getSVGElement());
        String ref = genome(2000);
        BaseRegion view = new BaseRegion(100, 130);
        BaseSeqViewModel refModel = BaseSeqViewModel.fromStr(ref.substring(89, 140), 90);
        Map<Integer, List<SvgRender.BoxBorder>> borders = new HashMap<>();
        borders.put(115, List.of(new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLACK), new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLACK)));
        borders.put(112, List.of(new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLUE)));
        borders.put(118, List.of(new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLUE)));
        borders.put(104, List.of(new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, ColorUtil.DARK_BLUE)));
        borders.put(126, List.of(new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, ColorUtil.DARK_BLUE)));
        out("svg.ref", SvgRender.renderBaseSeq(12.0, view, refModel, false, borders, null).getSVGElement());

        // a read: 3S8M2I6M3D9M2S starting (unclipped) at 99, bases from the reference with mismatches
        String cigar1 = "3S8M2I6M3D9M2S";
        StringBuilder rb = new StringBuilder();
        rb.append("TTT").append(ref, 101, 109).append("GG").append(ref, 109, 115).append(ref, 118, 127).append("NA");
        rb.setCharAt(6, rb.charAt(6) == 'A' ? 'C' : 'A');
        rb.setCharAt(20, rb.charAt(20) == 'G' ? 'T' : 'G');
        byte[] quals1 = new byte[rb.length()];
        for(int i = 0; i < quals1.length; i++) quals1[i] = (byte) (5 + (i * 7) % 36);
        List<CigarElement> cigarEls1 = TextCigarCodec.decode(cigar1).getCigarElements();
        out("read1.bases", rb);
        BaseSeqViewModel read1 = BaseSeqViewModel.create(99, cigarEls1, rb.toString().getBytes(), quals1, true);
        out("svg.read1", SvgRender.renderBaseSeq(12.0, view, read1, true, borders, refModel).getSVGElement());

        // a fragment: mates 20M at 96 and 5S18M at (unclipped) 108, the consensus 32M2S... formed by hand
        String first = ref.substring(95, 115);
        String second = "CCCCC" + ref.substring(112, 130);
        byte[] q30 = new byte[40];
        Arrays.fill(q30, (byte) 30);
        BaseSeqViewModel firstModel = BaseSeqViewModel.create(96, TextCigarCodec.decode("20M").getCigarElements(), first.getBytes(), q30, false);
        BaseSeqViewModel secondModel = BaseSeqViewModel.create(108, TextCigarCodec.decode("5S18M").getCigarElements(), second.getBytes(), q30, true);
        String consensus = ref.substring(95, 130);
        byte[] q37 = new byte[40];
        Arrays.fill(q37, (byte) 37);
        BaseSeqViewModel fragment = BaseSeqViewModel.create(96, TextCigarCodec.decode("35M").getCigarElements(), consensus.getBytes(), q37, false, firstModel, secondModel);
        out("svg.fragment", SvgRender.renderBaseSeq(12.0, view, fragment, true, borders, refModel).getSVGElement());
        out("svg.first", SvgRender.renderBaseSeq(12.0, view, firstModel, true, borders, refModel).getSVGElement());
        out("svg.second", SvgRender.renderBaseSeq(12.0, view, secondModel, true, borders, refModel).getSVGElement());

        // gene tracks from hand-made models
        List<GeneRegionViewModel> gm = new ArrayList<>();
        gm.add(new GeneRegionViewModel.IntronicRegionViewModel(new BaseRegion(90, 104), 0));
        gm.add(new GeneRegionViewModel.NonCodingExonicRegionViewModel(new BaseRegion(105, 110), 0));
        gm.add(new GeneRegionViewModel.AminoAcidViewModel(new BaseRegion(111, 113), 1, 'M', 'M', 0, true));
        gm.add(new GeneRegionViewModel.AminoAcidViewModel(new BaseRegion(114, 116), 2, 'W', 'K', 2, false));
        gm.add(new GeneRegionViewModel.DelViewModel(new BaseRegion(117, 119)));
        gm.add(new GeneRegionViewModel.AminoAcidViewModel(new BaseRegion(120, 122), 3, 'Q', '*', 0, false));
        gm.add(new GeneRegionViewModel.AminoAcidViewModel(new BaseRegion(123, 125), 3, 'I', 'I', 0, false));
        for(boolean pos : new boolean[] { true, false })
        {
            SvgRender.RenderedGeneData d = SvgRender.renderGeneData(12.0, new BaseRegion(100, 124), pos, gm, gm);
            out("svg.gene." + pos + ".ref", d.refSvgCanvas().getSVGElement());
            out("svg.gene." + pos + ".alt", d.altSvgCanvas().getSVGElement());
        }

        // HGVS parsing
        String[] hgvs = { "p.Gly12Val", "p.Gly12*", "p.Gly12=", "p.Gly12fs", "p.Met1?", "p.Gly12del", "p.Gly12_Ala14del", "p.Gly12dup",
                "p.Gly12_Ala14dup", "p.Gly12_Ala13insSerThr", "p.Gly12_Ala13insSer*", "p.Gly12_Ala14delinsVal", "p.Gly12delinsValSer",
                "p.Gly12_Ala14delinsValSerGly*", "p.Gly12_Ala14delinsValSerGlyThrPro", "p.Ter100Glnext*?", "p.Gly100ext*?", "p.Ter100Terext*?",
                "p.Gly12insSer", "p.Gly12Ter", "p.Gly12_Ala14delinsXyz", "p.12Gly", "c.12A>T", "p.", "p", "", "p.Gly12", "p.Gly12_Ala12del",
                "p.Gly12_Ter14del", "p.Ter12_Ala14del", "p.Gly12_Ala14fs", "p.Gly12_Ala14Val", "p.Gly12_Ala14insSer", "p.Gly12de", "p.Gly12delin",
                "p.Gly12delinsTer", "p.Gly12delinsVal*Ser", "p.Gly12_Ala13=", "p.Gly12_Ala13?", "p.Gly12_Ala14ext*?", "p.Gly12delins", "p.Gly12Xyz",
                "p.Ter100fs", "p.Ter100=", "p.Gly12_Ala14delinsValSer" };
        for(String h : hgvs)
        {
            String result;
            try { result = "ok " + events(AminoAcidEvent.parse(h)); }
            catch(Exception ex) { result = "err"; }
            out("hgvs." + h, result);
        }

        // gene regions and amino acid view models over a synthetic genome
        File fasta = new File(args[0]);
        try(PrintWriter w = new PrintWriter(fasta))
        {
            w.println(">chr1");
            for(int i = 0; i < ref.length(); i += 60) w.println(ref.substring(i, Math.min(ref.length(), i + 60)));
        }
        try(PrintWriter w = new PrintWriter(new File(args[0] + ".fai"))) { w.println("chr1\t2000\t6\t60\t61"); }
        try(PrintWriter w = new PrintWriter(new File(args[0] + ".dict"))) { w.println("@HD\tVN:1.6"); w.println("@SQ\tSN:chr1\tLN:2000"); }
        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(args[0]);
        String aminoAcids = "MACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVW";
        TranscriptAminoAcids taa = new TranscriptAminoAcids("ENSG0001", "GENE1", "ENST0001", true, aminoAcids);
        TranscriptAminoAcids shortTaa = new TranscriptAminoAcids("ENSG0001", "GENE1", "ENST0001", true, aminoAcids.substring(0, 30));
        Map<String, TranscriptData> transcripts = new LinkedHashMap<>();
        transcripts.put("pos", transcript((byte) 1, 1110, 1280));
        transcripts.put("neg", transcript((byte) -1, 1110, 1280));
        transcripts.put("posfull", transcript((byte) 1, 1100, 1300));
        transcripts.put("neg2", transcript((byte) -1, 1161, 1245));
        String[][] variants = { { "1180", "A", "C" }, { "1180", "ACGT", "A" }, { "1180", "A", "ACC" }, { "1128", "AC", "A" }, { "1205", "T", "TGGGG" }, { "1250", "ACGTACG", "A" } };
        String[] aaHgvs = { "p.Gly8Val", "p.Gly8*", "p.Gly8=", "p.Gly8fs", "p.Met1?", "p.Gly8del", "p.Gly6_Ala9del", "p.Gly8dup", "p.Gly8_Ala9insSerThr",
                "p.Gly7_Ala9delinsVal", "p.Gly8delinsValSer", "p.Gly7_Ala9delinsValSerGly*", "p.Ter38Glnext*?", "p.Gly20fs", "p.Gly30del", "p.Gly2fs",
                "p.Gly36_Ala39del", "p.Gly1_Ala3delinsValSerGlyThrPro" };
        for(Map.Entry<String, TranscriptData> entry : transcripts.entrySet())
        {
            TranscriptData t = entry.getValue();
            String name = entry.getKey();
            List<GeneRegionViewModel> regions = AminoAcidUtil.getGeneRegions(t);
            out("regions." + name, models(regions));
            out("refmodels." + name, models(AminoAcidUtil.getRefGeneRegionViewModels(t, taa, regions)));
            out("refmodels.short." + name, models(AminoAcidUtil.getRefGeneRegionViewModels(t, shortTaa, regions)));
            StringJoiner labels = new StringJoiner("|");
            for(int p : new int[] { 1000, 1100, 1115, 1130, 1131, 1159, 1160, 1210, 1220, 1240, 1300, 1301, 1400 })
                labels.add(AminoAcidUtil.getGeneRegionLabel(t, p));
            out("labels." + name, labels);
            for(String[] v : variants)
            {
                SimpleVariant sv = new SimpleVariant("chr1", Integer.parseInt(v[0]), v[1], v[2]);
                VariantReadContext rc = new VariantReadContext(sv, 0, 0, new byte[10], new byte[10], new ArrayList<>(), 1, 2, 3, null, null, new ArrayList<>(), 0, 0);
                SageVariant sageVariant = new SageVariant(new Candidate(VariantTier.PANEL, rc), new ArrayList<>(), new ArrayList<>());
                for(String h : aaHgvs)
                {
                    for(BaseRegion renderRegion : new BaseRegion[] { new BaseRegion(1105, 1281), new BaseRegion(1170, 1200) })
                    {
                        String key = "altmodels." + name + "." + v[0] + v[1] + ">" + v[2] + "." + h + "." + renderRegion.start();
                        try
                        {
                            List<GeneRegionViewModel> alt = AminoAcidUtil.getAltGeneRegionViewModels(t, taa, regions, AminoAcidEvent.parse(h), renderRegion, refGenome, sageVariant);
                            out(key, "ok " + models(alt));
                            if(renderRegion.start() == 1105 && (h.equals("p.Gly8fs") || h.equals("p.Gly8del") || h.equals("p.Gly8_Ala9insSerThr")))
                            {
                                SvgRender.RenderedGeneData d = SvgRender.renderGeneData(12.0, renderRegion, t.posStrand(), AminoAcidUtil.getRefGeneRegionViewModels(t, taa, regions), alt);
                                out("svg." + key, d.refSvgCanvas().getSVGElement() + d.altSvgCanvas().getSVGElement());
                            }
                        }
                        catch(Exception ex) { out(key, "err"); }
                    }
                }
            }
        }

        // read info table
        out("readinfo.paired", HtmlUtil.renderReadInfoTable("read<1>&\"x\"'", new ChrBaseRegion("chr1", 1000, 1099), new ChrBaseRegion("chr2", 1200, 1294),
                "100M", "5S95M", 295, "F1R2", 60, 3, "DUAL", "7;3;4", 59, 0, null).render());
        out("readinfo.single", HtmlUtil.renderReadInfoTable("r2", new ChrBaseRegion("chr1", 1000, 1049), null,
                "10S50M", "missing", 50, "", 13, 0, null, null, null, null, null).render());
        out("readinfo.nomc", HtmlUtil.renderReadInfoTable("r3", new ChrBaseRegion("chr1", 1000, 1049), new ChrBaseRegion("chr1", 900, 0),
                "50M", "missing", 150, "R1F2", 13, 2, "NONE", null, null, null, null).render());

        // HashSet<String> iteration order of filter names
        String[][] sets = { { "minTumorQual", "minTumorVAF", "minMQF", "maxGermlineRelQual" }, { "readStrandBias", "jitter", "minTumorQual", "minTumorVAF", "minMQF", "maxGermlineRelQual" },
                { "maxGermlineAltSupport", "minTumorQual", "minMQF", "maxGermlineRelQual", "maxGermlineVAF" }, { "dedupSnvMnv", "dedupIndel" } };
        for(String[] set : sets)
        {
            Set<String> hs = Arrays.stream(set).collect(Collectors.toSet());
            out("filters." + String.join(",", set), String.join(",", hs));
        }
    }
}
