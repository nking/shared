package algorithms.util;

import algorithms.sort.MultiArrayMergeSort;
import algorithms.util.LinkedListCostNode;

import java.math.BigInteger;
import java.util.logging.Logger;

import gnu.trove.map.TShortIntMap;
import gnu.trove.map.hash.TShortIntHashMap;
import gnu.trove.iterator.TShortIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.set.TShortSet;
import junit.framework.TestCase;

import java.util.Random;
import java.io.IOException;

/**
 * 
 * @author nichole
 */
public class FNVHashTest extends TestCase {

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

    private int oldHash(int[] params) {
        int fnv321aInit = 0x811c9dc5;
        int fnv32Prime = 0x01000193;
        int hash = fnv321aInit;
        for (int p : params) {
            hash = hash ^ p;
            hash = hash * fnv32Prime;
        }
        return hash;
    }

    public void testPairIntHashing() {
        /*
        there is a difference in the results of the Misra Gies edge coloring algorithm
        for use of PairInt.
        PairInt previously used the old fnv hash and I just changed it to a new one.
        the old hash produced better results.

        The classes which use the hash function in the Misra Gies Coloring.java are

            Set<PairInt> edges = new HashSet<PairInt>();
            and
            Map<PairInt, Integer> colorMap

        Here are the pairs of points for 2 different unit tests that produce different results:

        [junit] testEdgeColoring5
    [junit] (0,1)
    [junit] (0,2)
    [junit] (0,3)
    [junit] (0,1)
    [junit] (1,4)
    [junit] (1,5)
    [junit] (1,6)
    [junit] (0,2)
    [junit] (2,7)
    [junit] (2,8)
    [junit] (2,9)
    [junit] (0,3)
    [junit] (3,10)
    [junit] (3,11)
    [junit] (1,4)
    [junit] (1,5)
    [junit] (1,6)
    [junit] (2,7)
    [junit] (2,8)
    [junit] (2,9)
    [junit] (3,10)
    [junit] (3,11)
    [junit] testEdgeColoring5 k=4

    [junit] testEdgeColoring6
    [junit] (0,1)
    [junit] (0,2)
    [junit] (0,1)
    [junit] (1,3)
    [junit] (1,6)
    [junit] (0,2)
    [junit] (2,4)
    [junit] (2,6)
    [junit] (1,3)
    [junit] (3,5)
    [junit] (3,6)
    [junit] (2,4)
    [junit] (4,5)
    [junit] (4,6)
    [junit] (3,5)
    [junit] (4,5)
    [junit] (1,6)
    [junit] (2,6)
    [junit] (3,6)
    [junit] (4,6)
    [junit] testEdgeColoring6 k=5

    so will look at the differences in detail here, that is collisions, containskey
         */
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
        TIntIntMap freq2 = new TIntIntHashMap();

        final FNVHash fnv = new FNVHash();

        short limit = Short.MAX_VALUE;//32767

        Random rand = new Random();
        int r;
        short h;
        int h2;
        for (short i = 0; i < limit; ++i) {
            //r = rand.nextInt(limit);
            //h = fnv._hash(new short[]{r});

            h = fnv._hash(new short[]{i});
            if (freq.containsKey(h)) {
                freq.put(h, freq.get(h) + 1);
            } else {
                freq.put(h, 1);
            }

            h2 = oldHash(new int[]{i});
            if (freq2.containsKey(h2)) {
                freq2.put(h2, freq2.get(h2) + 1);
            } else {
                freq2.put(h2, 1);
            }
        }

        int[] xPoints = new int[freq.size()];
        int[] yPoints = new int[freq.size()];
        int[] xPoints2 = new int[freq2.size()];
        int[] yPoints2 = new int[freq2.size()];
        int[] xPolygon = null;
        int[] yPolygon = null;
        String label = "hash frequency";
        String label2 = "old hash frequency";

        TShortIntIterator iter = freq.iterator();
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
        plotter.writeFile("hash_freq");
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

        // key = hash, value = number of times seen
        TIntIntMap freq = new TIntIntHashMap();
        TIntIntMap freq2 = new TIntIntHashMap();

        final FNVHash fnv = new FNVHash();

        int limit = (1<<18)-1;

        Random rand = new Random();
        int r;
        int h;
        int h2;
        for (int i = 0; i < limit; ++i) {
            //r = rand.nextInt(limit);
            //h = fnv._hash(new short[]{r});

            h = fnv.hash(new int[]{i});
            if (freq.containsKey(h)) {
                freq.put(h, freq.get(h) + 1);
            } else {
                freq.put(h, 1);
            }

            h2 = oldHash(new int[]{i});
            if (freq2.containsKey(h2)) {
                freq2.put(h2, freq2.get(h2) + 1);
            } else {
                freq2.put(h2, 1);
            }
        }

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
    }

}
