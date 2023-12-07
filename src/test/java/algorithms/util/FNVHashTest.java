package algorithms.util;

import algorithms.misc.MiscMath0;
import algorithms.sort.MultiArrayMergeSort;
import algorithms.util.LinkedListCostNode;

import java.math.BigInteger;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

import gnu.trove.list.TFloatList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TShortIntMap;
import gnu.trove.map.hash.TShortIntHashMap;
import gnu.trove.iterator.TShortIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.set.TShortSet;
import junit.framework.TestCase;
import org.apache.datasketches.cpc.CpcSketch;
import org.apache.datasketches.frequencies.ErrorType;
import org.apache.datasketches.frequencies.ItemsSketch;
import org.apache.datasketches.kll.KllFloatsSketch;
import org.junit.Assert;

import java.util.Random;
import java.io.IOException;

/**
 * 
 * @author nichole
 */
public class FNVHashTest extends TestCase {


   /* see:
   https://datasketches.apache.org/docs/Community/Research.html
   https://datasketches.apache.org/docs/HLL/HLL.html
   https://datasketches.apache.org/docs/Sampling/ReservoirSampling.html
   */

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());
    
    public FNVHashTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testWithShort() throws IOException {
        log.info("testWithShort");

        // test hash for x bits
        // n = number of elements to hash
        //
        // expect 1st collision at sqrt(2)*sqrt(2^x)
        // expect all elements to collide at 2*(2^x)
        //   for 14 bits, 1st collision = 181, all elements collide = 32766
        //   for 15 bits, 1st collision = 256, all elements collide = 65534
        //   for 31 bits, 1st collision = 65536, all elements collide = 8.6e9

        // key = hash, value = number of times seen
        TShortIntMap freq = new TShortIntHashMap();
        PolygonAndPointPlotter plotter = null;//new PolygonAndPointPlotter();

        int limit = (1<<15)-1;
        // shorter for runTests. comment this out for full testing and see plotter above
        limit = (1 << 8) - 1;
        long n = 2*limit;
        final int lgK = 10;
        CpcSketch cpcSketch;
        KllFloatsSketch kLLSketch;
        ItemsSketch<Short> iSketch;

        Random rand = new Random();
        int r;
        short h;
        int key;
        String label;
        for (int type = 0; type < 2; ++type) {
            cpcSketch = new CpcSketch(lgK);
            kLLSketch = KllFloatsSketch.newHeapInstance();
            iSketch = new ItemsSketch<Short>(limit + 1);// size has to be a power of 2
            if (type == 0) {
                label = "old hash";
            } else {
                label = "new hash";
            }
            System.out.println(label);
            for (int i = 0; i < n; ++i) {
                r = rand.nextInt(limit);
                key = r;
                //key = i % limit;

                if (type == 0) {
                    h = FNVHash._oldhash(new short[]{(short)key});
                } else {
                    h = FNVHash._hash(new short[]{(short)key});
                }

                if (freq.containsKey(h)) {
                    freq.put(h, freq.get(h) + 1);
                } else {
                    freq.put(h, 1);
                }

                cpcSketch.update(h);
                kLLSketch.update(h);
                iSketch.update(h);
            }


            System.out.println(cpcSketch.toString());
            System.out.printf("Distinct count estimate: %f, expect=%d\n", cpcSketch.getEstimate(), limit);
            System.out.println("Distinct count lower bound 95% confidence: " + cpcSketch.getLowerBound(2));
            System.out.println("Distinct count upper bound 95% confidence: " + cpcSketch.getUpperBound(2));

            // hash min, max and median values
            float[] qs = kLLSketch.getQuantiles(new double[]{0.25, 0.5, 0.75, 1});
            System.out.printf("min=%f, max=%f, quantiles=%s\n", kLLSketch.getMinItem(), kLLSketch.getMaxItem(),
                    FormatArray.toString(qs, "%.3e"));

            // largest numbers of collisions
            //System.out.println("summary: " + iSketch.toString());
            //System.out.flush();
            ItemsSketch.Row<Short>[] items = iSketch.getFrequentItems(ErrorType.NO_FALSE_POSITIVES);
            System.out.println("number of Largest collisions (Frequent items): " + items.length);
            System.out.println(ItemsSketch.Row.getRowHeader());
            for (ItemsSketch.Row<Short> row : items) {
                System.out.println("First item: " + row.toString());
                break;
            }

            int[] xPoints = new int[freq.size()];
            int[] yPoints = new int[freq.size()];
            int[] xPolygon = null;
            int[] yPolygon = null;

            TShortIntIterator iter = freq.iterator();
            int i = 0;
            while (iter.hasNext()) {
                iter.advance();
                xPoints[i] = iter.key();
                yPoints[i] = iter.value();
                ++i;
            }
            MultiArrayMergeSort.sortByIncr(xPoints, yPoints);

            // plot the hash vs frequency.
            if (plotter != null) {
                plotter.addPlot(xPoints, yPoints, xPolygon, yPolygon, label);
            }
            // FreedmanDiaconis: binWidth = 2*IQR * n^(−1/3)
            long nn = kLLSketch.getN();
            float iqr = kLLSketch.getQuantile(0.75) - kLLSketch.getQuantile(0.5);
            double binWidth = 2. * iqr * Math.pow(nn, -0.33333);
            double dwd = binWidth/nn;
            binWidth = 2. * iqr * Math.pow(limit, -0.33333);
            dwd = binWidth/limit;

            // plot rank vs quantile
            float dw = 0.01f;
            int nw = (int)(1.f/dw);
            TFloatList xd = new TFloatArrayList();
            TFloatList yd = new TFloatArrayList();
            for (i = 0; i < nw; ++i) {
                xd.add(i*dw);
                yd.add(kLLSketch.getQuantile(xd.get(i)));
            }
            /*double[] pmf = kLLSketch.getPMF(xd);
            for (i = 0; i < nw; ++i) {
                yd[i] = (float)pmf[i];
            }*/
            if (plotter != null) {
                plotter.addPlot(xd.toArray(), yd.toArray(), xPolygon, yPolygon, label + " CDF");
            }
            // create splitpoints and PMF
            dw = (float)binWidth;
            nw = (int)(limit/dw);
            xd.clear();
            yd.clear();
            xd.add(1.f);
            for (i = 0; i < nw; ++i) {
                xd.add(i*dw);
            }
            double[] pmf = kLLSketch.getPMF(xd.toArray());
            double sum = 0;
            for (i = 0; i < nw; ++i) {
                yd.add((float)pmf[i]);
                sum += yd.get(i);
            }
            if (plotter != null) {
                plotter.addPlot(xd.toArray(), yd.toArray(), xPolygon, yPolygon, label + " PMF.  sum=" + sum);
            }
        }
        if (plotter != null) {
            plotter.writeFile("hash_freq");
        }

        /*
        int limit = (1<<15)-1;
        long n = 2*limit;

        old hash
        number of Largest collisions (Frequent items): 18950
                   Est          UB          LB Item
        First item:             24          24          24 16138
        new hash
        number of Largest collisions (Frequent items): 18939
                   Est          UB          LB Item
        First item:             26          26          26 15908

         */
    }

    public void testWithInt() throws IOException {
        log.info("testWithInt");

        // test hash for x bits
        // n = number of elements to hash
        //
        // expect 1st collision at sqrt(2)*sqrt(2^x)
        // expect all elements to collide at 2*(2^x)
        //   for 14 bits, 1st collision = 181, all elements collide = 32766
        //   for 15 bits, 1st collision = 256, all elements collide = 65534
        //   for 31 bits, 1st collision = 65536, all elements collide = 8.6e9

        int limit = (1<<24)-1;//(1<<30)-1;
        // shorter for runTests
        limit = (1 << 8) - 1;
        long n = 2*limit;
        final int lgK = 10;
        //CpcSketch cpcSketch;
        //KllFloatsSketch kLLSketch;
        ItemsSketch<Integer> iSketch;

        Random rand = new Random();
        int r;
        int h;
        int key;
        for (int type = 0; type < 2; ++type) {
            //cpcSketch = new CpcSketch(lgK);
            //kLLSketch = KllFloatsSketch.newHeapInstance();
            iSketch = new ItemsSketch<Integer>(limit+1);// size has to be a power of 2
            if (type == 0) {
                System.out.println("old hash");
            } else {
                System.out.println("new hash");
            }
            for (int i = 0; i < n; ++i) {
                r = rand.nextInt(limit);
                key = r;
                //key = i;

                if (type == 0) {
                    h = FNVHash.hash32a(new int[]{key});
                } else {
                    h = FNVHash.hash31a(new int[]{key});
                }

                //cpcSketch.update(h);
                //kLLSketch.update(h);
                iSketch.update(h);
            }

            /*
            System.out.println(cpcSketch.toString());
            System.out.printf("Distinct count estimate: %f, expect=%d\n", cpcSketch.getEstimate(), limit);
            System.out.println("Distinct count lower bound 95% confidence: " + cpcSketch.getLowerBound(2));
            System.out.println("Distinct count upper bound 95% confidence: " + cpcSketch.getUpperBound(2));

            // hash min, max and median values
            System.out.printf("min=%f, max=%f, median=%f\n", kLLSketch.getMinItem(), kLLSketch.getMaxItem(),
                    kLLSketch.getQuantile(0.5));
            */

            // largest numbers of collisions
            //System.out.println("summary: " + iSketch.toString());
            //System.out.flush();
            ItemsSketch.Row<Integer>[] items = iSketch.getFrequentItems(ErrorType.NO_FALSE_POSITIVES);
            System.out.println("number of Largest collisions (Frequent items): " + items.length);
            System.out.println(ItemsSketch.Row.getRowHeader());
            for (ItemsSketch.Row<Integer> row: items) {
                System.out.println("First item: " + row.toString());
                break;
            }
        }
        /*
        with int limit = (1<<22)-1;//(1<<30)-1;
        long n = 2*limit;

         old hash
        Largest collisions (Frequent items): 876563
                   Est          UB          LB Item
        First item:             12          12          10 -2146108584

        new hash
        Largest collisions (Frequent items): 875945
                   Est          UB          LB Item
        First item:             14          14          12 434781249

        -----
        with int limit = (1<<24)-1;//(1<<30)-1;
        long n = 2*limit;

        old hash
        number of Largest collisions (Frequent items): 3506275
                   Est          UB          LB Item
        First item:             13          13          11 780217474
        new hash
        number of Largest collisions (Frequent items): 3497277
                   Est          UB          LB Item
        First item:             14          14          12 317334821
         */
        /*
        int[] xPoints = new int[freq.size()];
        int[] yPoints = new int[freq.size()];
        int[] xPoints2 = new int[freq2.size()];
        int[] yPoints2 = new int[freq2.size()];
        int[] xPolygon = null;
        int[] yPolygon = null;
        String label = "hash frequency int";
        String label2 = "old hash frequency int";

        TIntIntIterator iter = freq.iterator();
        int i = 0;
        while (iter.hasNext()) {
            iter.advance();
            xPoints[i] = iter.key();
            yPoints[i] = iter.value();
            ++i;
        }
        MultiArrayMergeSort.sortByIncr(xPoints, yPoints);
        TIntIntIterator iter2 = freq2.iterator();
        i = 0;
        while (iter2.hasNext()) {
            iter2.advance();
            xPoints2[i] = iter2.key();
            yPoints2[i] = iter2.value();
            ++i;
        }
        MultiArrayMergeSort.sortByIncr(xPoints2, yPoints2);
        // plot the hash vs frequency.
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(xPoints, yPoints, xPolygon, yPolygon, label);
        plotter.addPlot(xPoints2, yPoints2, xPolygon, yPolygon, label2);
        plotter.writeFile("hash_freq_int");
        */
    }

}
