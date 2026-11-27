package com.hartwig.hmftools.bamtools.metrics;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

// Synthetic cases through Java's CoverageMetrics statistics, coverage buckets and OffTargetAnalyser.
// Each case writes its input and Java's output on one line, for the Rust test to replay.
public class BtHarness
{
    public static void main(String[] args) throws Exception
    {
        PrintWriter out = new PrintWriter(args[0]);
        Random rnd = new Random(42);

        // statistics: random depth distributions
        for(int c = 0; c < 400; ++c)
        {
            int maxCov = 1 + rnd.nextInt(c < 300 ? 300 : 12000);
            CoverageMetrics m = new CoverageMetrics(maxCov);
            int shape = c % 4;
            for(int i = 0; i <= maxCov; ++i)
            {
                long f;
                if(shape == 0) f = rnd.nextInt(3) == 0 ? rnd.nextInt(1000) : 0;
                else if(shape == 1) f = (long)(1e6 * Math.exp(-Math.pow((i - maxCov / 3.0) / (maxCov / 10.0 + 1), 2)));
                else if(shape == 2) f = (i < 3 && rnd.nextBoolean()) ? rnd.nextInt(5) : 0;
                else f = rnd.nextInt(100000);
                m.CoverageFrequency[i] = f;
            }
            for(int t = 0; t < m.FilterTypeCounts.length; ++t)
                m.FilterTypeCounts[t] = rnd.nextInt(4) == 0 ? 0 : (long)rnd.nextInt(Integer.MAX_VALUE) * rnd.nextInt(1000);
            boolean excludeZero = rnd.nextBoolean();
            long zero = rnd.nextInt(1000000), cov = rnd.nextInt(100000000);
            m.addCoverageBases(cov, zero);
            StringBuilder sb = new StringBuilder("STATS\t" + maxCov + "\t" + excludeZero + "\t" + zero + "\t" + cov + "\t");
            for(int i = 0; i <= maxCov; ++i) sb.append(i > 0 ? "," : "").append(m.CoverageFrequency[i]);
            sb.append("\t");
            for(int t = 0; t < m.FilterTypeCounts.length; ++t) sb.append(t > 0 ? "," : "").append(m.FilterTypeCounts[t]);
            m.finalise(excludeZero);
            Statistics s = m.statistics();
            sb.append("\t").append(Double.doubleToLongBits(s.Mean)).append("\t").append(Double.doubleToLongBits(s.Median))
                    .append("\t").append(s.MedianAbsoluteDeviation).append("\t").append(Double.doubleToLongBits(s.StandardDeviation));
            sb.append("\t");
            for(FilterType t : FilterType.values()) sb.append(t.ordinal() > 0 ? "," : "").append(Double.doubleToLongBits(m.calcFilteredPercentage(t)));
            sb.append("\t");
            for(int k = 0; k < MetricsWriter.COVERAGE_LEVELS.size(); ++k)
                sb.append(k > 0 ? "," : "").append(Double.doubleToLongBits(m.calcCoverageFrequency(MetricsWriter.COVERAGE_LEVELS.get(k), m.genomeTerritory())));
            out.println(sb);
        }

        // coverage buckets
        for(int n : new int[] { 3001, 5000, 10001, 20001 })
        {
            int[] buckets = CoverageMetrics.createCoverageBuckets(n);
            StringBuilder sb = new StringBuilder("BUCKETS\t" + n + "\t");
            for(int i = 0; i < buckets.length; ++i) sb.append(i > 0 ? "," : "").append(buckets[i]);
            sb.append("\t");
            for(int i = 1; i < n; ++i) sb.append(i > 1 ? "," : "").append(CoverageMetrics.getNextCoverageBucket(i, buckets));
            out.println(sb);
        }

        // off-target pile-ups: sorted fragment starts, random lengths, the same read now and then
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 100000000));
        for(int c = 0; c < 300; ++c)
        {
            int threshold = rnd.nextInt(6);
            OffTargetAnalyser a = new OffTargetAnalyser(threshold);
            int n = 1 + rnd.nextInt(60);
            int pos = 1000;
            StringBuilder in = new StringBuilder();
            List<String> names = new ArrayList<>();
            for(int i = 0; i < n; ++i)
            {
                pos += rnd.nextInt(rnd.nextInt(4) == 0 ? 400 : 60);
                int len = 20 + rnd.nextInt(300);
                String name = (rnd.nextInt(10) == 0 && !names.isEmpty()) ? names.get(rnd.nextInt(names.size())) : "r" + i;
                names.add(name);
                SAMRecord r = new SAMRecord(header);
                r.setReadName(name);
                r.setReferenceName("chr1");
                r.setAlignmentStart(pos);
                r.setCigarString(len + "M");
                r.setReadPairedFlag(false);
                a.addRead(r, true, 0);
                in.append(i > 0 ? ";" : "").append(name).append(",").append(pos).append(",").append(pos + len - 1);
            }
            List<String> sites = a.enrichedFragmentSites();
            Map<Integer,Integer> counts = a.fragmentOverlapCounts();
            List<Integer> keys = new ArrayList<>(counts.keySet());
            Collections.sort(keys);
            StringBuilder sb = new StringBuilder("OFFTARGET\t" + threshold + "\t" + in + "\t");
            for(int k = 0; k < keys.size(); ++k) sb.append(k > 0 ? "," : "").append(keys.get(k)).append("=").append(counts.get(keys.get(k)));
            sb.append("\t").append(String.join("|", sites));
            out.println(sb);
        }
        out.close();
    }
}
