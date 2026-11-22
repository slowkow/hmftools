import java.awt.Color;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.*;
import org.jfree.svg.SVGGraphics2D;
import org.apache.commons.lang3.tuple.Pair;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.vis.*;
import com.hartwig.hmftools.common.vis.SvgRender.ChrLabel;
import com.hartwig.hmftools.common.vis.SvgUtil.Alignment;
import static j2html.TagCreator.*;

// Prints expected values for the Rust unit tests of the copied vis helpers, from synthetic input only.
public class VisProbe {
    static void out(String name, String value) { System.out.println(name + "\t" + value); }
    static List<CigarElement> cigar(String s) { return new ArrayList<>(htsjdk.samtools.TextCigarCodec.decode(s).getCigarElements()); }

    public static void main(String[] a) {
        // CssBuilder / CssSize
        CssBuilder border = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder firstRow = CssBuilder.EMPTY.borderLeft(CssSize.px(1.0), "solid", Color.BLACK).borderRight(CssSize.px(1.0), "solid", Color.BLACK)
                .borderTop(CssSize.px(2.0), "solid", Color.BLACK).borderBottom(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        out("css.empty", "[" + CssBuilder.EMPTY.toString() + "]");
        out("css.table", border.marginRight(CssSize.px(10.0)).toString());
        out("css.cell", CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(border).textAlign("right").toString());
        out("css.left_header", CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).writingMode("vertical-rl").transform("rotate(180deg)").fontWeight("bold").textAlign("center").merge(firstRow).toString());
        out("css.merge_wins", CssBuilder.EMPTY.textAlign("left").height(CssSize.em(1)).merge(CssBuilder.EMPTY.textAlign("right")).toString());
        out("css.div", CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO).display("none").toString());
        out("css.base_font", HtmlUtil.BASE_FONT_STYLE.toString());
        out("css.sizes", CssSize.px(0.00004) + " " + CssSize.px(0.00005) + " " + CssSize.em(1.23456) + " " + CssSize.px(-2.5) + " " + CssSize.ZERO);
        out("css.shadow", CssBuilder.EMPTY.leftBoxShadow(CssSize.px(2), Color.RED).rightBoxShadow(CssSize.px(3), Color.BLUE).fontSizePct(0.755).toString());

        // j2html
        out("html.escape", td("a<b>&\"c'd").render());
        out("html.attrs", td("x").attr("rowspan", 3).withStyle("color:#000000;").render());
        out("html.table", HtmlUtil.styledTable(List.of(tr(td(""), td(rawHtml("<svg/>")))), CssBuilder.EMPTY.borderSpacing(CssSize.ZERO)).render());
        String js = HtmlUtil.getJavascript().render();
        out("html.js", js.length() + " " + js.hashCode() + " " + js.substring(0, 40).replace("\n", "\\n") + " " + js.substring(js.length() - 20).replace("\n", "\\n"));

        // read info table
        List<Pair<String, String>> extra = new ArrayList<>();
        extra.add(Pair.of("Mismatch Info:", "1:2 <x>"));
        extra.add(Pair.of("Trim count:", "0,3"));
        out("info.full", HtmlUtil.renderReadInfoTable("read<1>", new ChrBaseRegion("chrA", 100, 250), new ChrBaseRegion("chrB", 400, 0), "10S140M", null, 345, null, 60, null, null, null, null, null, extra).render());
        out("info.unmapped", HtmlUtil.renderReadInfoTable("read2", new ChrBaseRegion("chrA", 100, 250), null, "150M", "100M", null, "FR", null, 3, "DUAL", "5;2;1", 20, 4, null).render());

        // string bounds
        for(String s : new String[] { "9", "1,234,567", ">> chr7:1234567", "chrUn_1:99 <<", "INS", "A", "128" })
            for(java.awt.Font f : new java.awt.Font[] { SvgRender.COORD_FONT, new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 26), new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 18) }) {
                java.awt.geom.Rectangle2D r = SvgUtil.getStringBounds(f, s);
                out("bounds." + f.getSize() + "." + s, r.getX() + " " + r.getWidth() + " " + r.getHeight() + " " + r.getCenterX() + " " + SvgUtil.getFontDescent(f));
            }

        // base sequence view models
        byte[] bases = "ACGTTGCANNACGTACGTAC".getBytes();
        byte[] quals = new byte[bases.length];
        for(int i = 0; i < quals.length; i++) quals[i] = (byte) (i * 2 + 1);
        quals[5] = (byte) 200; // negative as a Java byte
        BaseSeqViewModel ref = BaseSeqViewModel.fromStr("ACGTAGCATTACGGACGAACGTTT", 100);
        BaseSeqViewModel read = BaseSeqViewModel.create(101, cigar("2S4M2I3M3D5M1H4S"), bases, quals, true);
        out("model.read", read.FirstBasePos + " " + read.LastBasePos + " " + read.LeftIsForwardStrand + " " + read.RightIsForwardStrand + " " + read.mismatchCount(ref) + " " + read.getBase(101).isSoftClip() + " " + read.getBase(106).rightInsertCount() + " " + read.getBase(110).isDel() + " " + read.getBase(118).isMissing() + " " + read.getBase(99).isMissing() + " " + read.getBase(104).baseQ());
        read.clearSoftClips();
        out("model.read_cleared", read.mismatchCount(ref) + " " + read.getBase(101).isSoftClip());
        BaseSeqViewModel asm = BaseSeqViewModel.fromStringWithCigar("ACGTAGCTTTACGG", cigar("3S5M1D2I4M2N"), 98);
        out("model.asm", asm.FirstBasePos + " " + asm.LastBasePos + " " + asm.hasOrientation() + " " + asm.mismatchCount(ref) + " " + asm.getBase(98).isSoftClip() + " " + asm.getBase(106).isDel() + " " + asm.getBase(105).rightInsertCount());
        BaseSeqViewModel empty = BaseSeqViewModel.fromStr("", 5);
        out("model.empty", empty.FirstBasePos + " " + empty.LastBasePos + " " + read.mismatchCount(empty));
        try { BaseSeqViewModel.create(1, cigar("30M"), bases, quals, false); out("model.overrun", "none"); } catch(Exception e) { out("model.overrun", e.getClass().getName()); }
        try { BaseSeqViewModel.fromStringWithCigar("ACGT", cigar("5M"), 1); out("model.str_overrun", "none"); } catch(Exception e) { out("model.str_overrun", e.getClass().getName()); }
        try { BaseSeqViewModel m = BaseSeqViewModel.create(1, cigar("3I2M"), "ACGTA".getBytes(), new byte[]{1,2,3,4,5}, false); out("model.lead_insert", m.FirstBasePos + " " + m.LastBasePos + " " + m.getBase(1).charBase()); } catch(Exception e) { out("model.lead_insert", e.getClass().getName()); }

        // render base sequences on one canvas at two offsets, as the assembly visualiser does
        BaseSeqViewModel read2 = BaseSeqViewModel.create(101, cigar("2S4M2I3M3D5M1H4S"), bases, quals, true);
        SVGGraphics2D canvas = new SVGGraphics2D(12.0 * 40, 12.0);
        AffineTransform init = canvas.getTransform();
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(0.0, 0.0), 12.0, new BaseRegion(100, 112), read2, true, new HashMap<>(), ref, true, false, false);
        canvas.setTransform(init);
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(15.0, 0.0), 12.0, new BaseRegion(108, 123), read2, true, new HashMap<>(), ref, false, true, true);
        out("svg.read", canvas.getSVGElement());

        canvas = new SVGGraphics2D(12.0 * 30, 12.0);
        init = canvas.getTransform();
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(0.0, 0.0), 12.0, new BaseRegion(100, 110), ref, false, new HashMap<>(), null);
        canvas.setTransform(init);
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(13.0, 0.0), 12.0, new BaseRegion(98, 112), asm, false, new HashMap<>(), ref);
        out("svg.ref_asm", canvas.getSVGElement());

        BaseSeqViewModel fwd = BaseSeqViewModel.create(100, cigar("5M"), "ACGTT".getBytes(), new byte[] { 37, 40, 0, 10, 20 }, false);
        canvas = new SVGGraphics2D(12.0 * 9, 12.0);
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(1.0, 0.0), 12.0, new BaseRegion(100, 104), fwd, true, new HashMap<>(), ref, true, true, false);
        out("svg.fwd", canvas.getSVGElement());
        canvas = new SVGGraphics2D(12.0 * 9, 12.0);
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(1.0, 0.0), 12.0, new BaseRegion(100, 104), fwd, true, new HashMap<>(), ref, true, true, true);
        out("svg.fwd_reversed_markers", canvas.getSVGElement());
        // a deletion that runs to the end of the region, and an insert before the region
        BaseSeqViewModel del = BaseSeqViewModel.create(100, cigar("2M2I1M6D2M"), "ACGGTAC".getBytes(), new byte[] { 30, 30, 30, 30, 30, 30, 30 }, false);
        canvas = new SVGGraphics2D(12.0 * 8, 12.0);
        SvgRender.renderBaseSeq(canvas, new Point2D.Double(0.0, 0.0), 12.0, new BaseRegion(102, 105), del, false, new HashMap<>(), ref, true, true, false);
        out("svg.del_end", canvas.getSVGElement());

        // coordinates: two regions on one canvas, the second reversed
        double scalingFactor = 12.0 / SvgRender.BASE_BOX_SIZE;
        double charWidth = SvgUtil.getStringBounds(SvgRender.COORD_FONT, "9").getWidth();
        double maxStringWidth = SvgUtil.getStringBounds(SvgRender.COORD_FONT, "1,000,012").getWidth();
        BaseRegion r1 = new BaseRegion(999990, 1000012);
        BaseRegion r2 = new BaseRegion(2031, 2050);
        int boxes = r1.baseLength() + 2 + r2.baseLength() + 2;
        canvas = new SVGGraphics2D(scalingFactor * SvgRender.BASE_BOX_SIZE * boxes, scalingFactor * (maxStringWidth + 2 * charWidth));
        init = canvas.getTransform();
        SvgRender.renderCoords(canvas, new Point2D.Double(0.0, 0.0), new Point2D.Double(scalingFactor * SvgRender.BASE_BOX_SIZE * (r1.baseLength() + 2), canvas.getHeight()), 12.0, r1, 1000012, 10, maxStringWidth, false);
        canvas.setTransform(init);
        SvgRender.renderCoords(canvas, new Point2D.Double(r1.baseLength() + 2, 0.0), new Point2D.Double(scalingFactor * SvgRender.BASE_BOX_SIZE * (r2.baseLength() + 2), canvas.getHeight()), 12.0, r2, 2031, 10, maxStringWidth, true);
        out("svg.coords", canvas.getSVGElement());

        // chromosome labels with an insert between
        List<ChrLabel> labels = new ArrayList<>();
        labels.add(new ChrLabel("chr7", 1500, new BaseRegion(0, 24), new BaseRegion(1001, 1500), false, List.of(Alignment.RIGHT)));
        labels.add(new ChrLabel(null, 0, new BaseRegion(25, 31), null, false, null));
        labels.add(new ChrLabel("12", 88000, new BaseRegion(32, 60), new BaseRegion(88000, 88400), true, List.of(Alignment.LEFT)));
        out("svg.chr_labels", SvgRender.renderChrLabels(12.0, labels).getSVGElement());
        labels = new ArrayList<>();
        labels.add(new ChrLabel("chrX", 20, new BaseRegion(0, 30), new BaseRegion(20, 300), true, List.of(Alignment.RIGHT)));
        labels.add(new ChrLabel("chrX", 900, new BaseRegion(31, 70), new BaseRegion(600, 900), false, List.of(Alignment.LEFT)));
        out("svg.chr_labels2", SvgRender.renderChrLabels(12.0, labels).getSVGElement());

        // PriorityQueue poll order with ties: elements are {key1, key2, tag}
        int[][] items = { {5,9,0}, {3,7,1}, {5,9,2}, {3,7,3}, {1,2,4}, {5,9,5}, {3,7,6}, {8,8,7}, {1,2,8}, {5,1,9}, {3,7,10}, {0,0,11}, {5,9,12} };
        PriorityQueue<int[]> queue = new PriorityQueue<>(Comparator.comparingInt((int[] x) -> x[0]).thenComparingInt((int[] x) -> x[1]));
        queue.addAll(Arrays.asList(items));
        StringBuilder sb = new StringBuilder();
        while(!queue.isEmpty()) sb.append(queue.poll()[2]).append(' ');
        out("pq.order", sb.toString().trim());

        // colours
        Color c = ColorUtil.lighten(new Color(151, 83, 9), 1.0 - 1.0 * 20 / 37);
        out("color.lighten", c.getRed() + "," + c.getGreen() + "," + c.getBlue() + " " + SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR.getRed() + "," + SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR.getGreen() + "," + SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR.getBlue());
    }
}
