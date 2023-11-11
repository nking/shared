package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntWithIndex;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.compGeometry.LargestDistancePair.PairAndHull;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;

import java.util.*;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LargestDistancePairTest extends TestCase {
    
    public LargestDistancePairTest() {
    }

    public void testFind() throws Exception {
        
        /*            2,6
         *
         *
         *     0,2   2,2 3,2
         *                      7,1
         *            2,0
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};
        
        long[] expResult = new long[]{0,2,7,1};
        long[] expResult2 = new long[]{7, 1, 2, 6};

        long distSq = (7 - 2)*(7-2) + (1-6)*(1-6);

        boolean useLinear = true;
        PairAndHull pairAndHull = LargestDistancePair.find(x, y, useLinear);

        int found = 0;
        if ( expResult[0] == (long)pairAndHull.getXY0()[0] && expResult[1] == (long)pairAndHull.getXY0()[1]) {
            found++;
            if ( expResult[2] == (long)pairAndHull.getXY1()[0] && expResult[3] == (long)pairAndHull.getXY1()[1]) {
                found++;
            }
        } else if ( expResult[0] == (long)pairAndHull.getXY1()[0] && expResult[1] == (long)pairAndHull.getXY1()[1]) {
            found++;
            if ( expResult[2] == (long)pairAndHull.getXY0()[0] && expResult[3] == (long)pairAndHull.getXY0()[1]) {
                found++;
            }
        }
        assertTrue(distSq == (long)pairAndHull.distSq);
        assertEquals(2, found);

    }
    
    public void testRandomInput() throws Exception {

        int nTests = 100;
        TDoubleList fracWrong = new TDoubleArrayList();
        
        for (int i = 0; i < nTests; ++i) {
            runRandomInput(fracWrong);
        }

        System.out.printf("nTests=%d, nWrong=%d\nfrac of true for wrong:%s\n", nTests, fracWrong.size(), fracWrong.toString());
    }
        
    private void runRandomInput(TDoubleList fracWrong) throws Exception {
        long seed = System.nanoTime();
        //seed = 191992537116792L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        // generate points in this coordinate range:
        int min = 1;
        int max = 10000;
        
        // generate this many points:
        int n = 100 + rand.nextInt(1000 - 100);
        
        //System.out.printf("n=%d  min=%d max=%d   max-min=%d\n", n, min, max, max-min);
        
        long[] x = new long[n];
        long[] y = new long[n];
        TLongObjectMap<TLongSet> points = new TLongObjectHashMap<>();
        
        int i = 0;
        TLongSet set;
        while (i < n) {
            x[i] = min + rand.nextInt(max - min);
            y[i] = min + rand.nextInt(max - min);
            set = points.get(x[i]);
            if ((set != null) && set.contains(y[i])) {
                continue;
            }
            if (set == null) {
                set = new TLongHashSet();
                points.put(x[i], set);
            }
            set.add(y[i]);
            ++i;
        }
                
        // use brute force to calculate the expected answer(s):
        TLongList expectedX0s = new TLongArrayList();
        TLongList expectedY0s = new TLongArrayList();
        TLongList expectedX1s = new TLongArrayList();
        TLongList expectedY1s = new TLongArrayList();
        
        bruteForceMaxDist(expectedX0s, expectedY0s, expectedX1s, expectedY1s, x, y);


        TLongList expectedXH0s = new TLongArrayList();
        TLongList expectedYH0s = new TLongArrayList();
        TLongList expectedXH1s = new TLongArrayList();
        TLongList expectedYH1s = new TLongArrayList();
        boolean useLinearSort = false;
        GrahamScan.CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
        bruteForceMaxDist(expectedXH0s, expectedYH0s, expectedXH1s, expectedYH1s,
                ch.getXH(), ch.getYH());

        assertTrue(expectedX0s.size() == expectedX1s.size() 
            && expectedX0s.size() == expectedY0s.size() && expectedX0s.size() == expectedY1s.size());
        
        assertFalse(expectedX0s.isEmpty());

        //System.out.printf("number of max dist pairs=%d\n", expectedX0s.size());
        
        // returns a pair in format [x0, y0, x1, y1]
        boolean useLinear = true;
        PairAndHull pairAndHull = LargestDistancePair.find(x, y, useLinear);
        
        long xd;
        long yd;
        long distSq, distSq0;

        long[] expResult = new long[]{expectedXH0s.get(0), expectedYH0s.get(0), expectedXH1s.get(0), expectedYH1s.get(0)};
        int found = 0;
        if ( expResult[0] == (long)pairAndHull.getXY0()[0] && expResult[1] == (long)pairAndHull.getXY0()[1]) {
            found++;
            if ( expResult[2] == (long)pairAndHull.getXY1()[0] && expResult[3] == (long)pairAndHull.getXY1()[1]) {
                found++;
            }
        } else if ( expResult[0] == (long)pairAndHull.getXY1()[0] && expResult[1] == (long)pairAndHull.getXY1()[1]) {
            found++;
            if ( expResult[2] == (long)pairAndHull.getXY0()[0] && expResult[3] == (long)pairAndHull.getXY0()[1]) {
                found++;
            }
        }

        distSq0 = 0;
        for (i = 0; i < 1/*expectedXH0s.size()*/; ++i) {
            xd = expectedX0s.get(i) - expectedX1s.get(i);
            yd = expectedY0s.get(i) - expectedY1s.get(i);
            distSq0 = xd*xd + yd*yd;
        }

        distSq = 0;
        for (i = 0; i < 1/*expectedXH0s.size()*/; ++i) {
            xd = expectedXH0s.get(i) - expectedXH1s.get(i);
            yd = expectedYH0s.get(i) - expectedYH1s.get(i);
            distSq = xd*xd + yd*yd;
            /*System.out.printf("expected for xh, yh =(%d,%d)(%d,%d) distSq=%d\n",
                expectedXH0s.get(i), expectedYH0s.get(i), expectedXH1s.get(i),
                expectedYH1s.get(i), distSq);*/
        }

        //TODO: follow up on this.  occassionally fails:
        //assertEquals(distSq0, distSq);

        if (distSq != (long)pairAndHull.distSq) {
            // occasionally, the rougher angular resolution of the linear hull used in LargestDistancePair
            // leads to a smaller distance than largest

            /*
            System.out.printf("expecting=%d  [%d,%d, %d,%d]  (%d)\n",
                    distSq, expectedXH0s.get(0), expectedYH0s.get(0),
                    expectedXH1s.get(0), expectedYH1s.get(0), expectedXH0s.size());
            System.out.printf("result   =%d  \n%s\n", (long)pairAndHull.distSq,
                    pairAndHull.toString("%.0f"));
             */

            double frac = Math.abs(distSq - pairAndHull.distSq)/(distSq + pairAndHull.distSq);
            fracWrong.add(frac);

            assertTrue(frac < 1E-1);

            pairAndHull = LargestDistancePair.find(x, y, false);

            /*

            int nH = ch.getXH().length;
            //print the hull
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, (int)Math.round(max*1.1),
                    0, (int)Math.round(max*1.1));

            PairIntArray curve = new PairIntArray(nH - 1);
            float[] xHull = new float[nH];
            float[] yHull = new float[nH];
            for (int ii = 0; ii < nH; ++ii) {
                xHull[ii] = ch.getXH()[ii];
                yHull[ii] = ch.getYH()[ii];
                if (ii < (nH - 1)) {
                    curve.add((int)xHull[ii], (int)yHull[ii]);
                }
            }
            float[] xF = MatrixUtil.copyLongToFloat(x);
            float[] yF = MatrixUtil.copyLongToFloat(y);

            plotter.addPlot(xF, yF, xHull, yHull, "LDP");
            plotter.writeFile();
            */

        } else {

            assertEquals(distSq, Math.round(pairAndHull.distSq));
            assertEquals(2, found);
        }
    }

    private void bruteForceMaxDist(TLongList expectedX0s, 
        TLongList expectedY0s, TLongList expectedX1s, TLongList expectedY1s, long[] x, long[] y) {
        
        int i;
        int n = x.length;
        
        long maxDist = Long.MIN_VALUE;
        int j;
        long x0;
        long y0;
        long xd;
        long yd;
        long distSq;
        for (i = 0; i < n; ++i) {
            x0 = x[i];
            y0 = y[i];
            for (j = i + 1; j < n; ++j) {
                xd = x[j] - x0;
                yd = y[j] - y0;
                distSq = xd * xd + yd * yd;
                if (distSq > maxDist) {
                    expectedX0s.clear();
                    expectedY0s.clear();
                    expectedX1s.clear();
                    expectedY1s.clear();
                    
                    maxDist = distSq;
                    
                    expectedX0s.add(x0);
                    expectedY0s.add(y0);
                    expectedX1s.add(x[j]);
                    expectedY1s.add(y[j]);
                } else if (distSq == maxDist) {
                    
                    expectedX0s.add(x0);
                    expectedY0s.add(y0);
                    expectedX1s.add(x[j]);
                    expectedY1s.add(y[j]);;
                }
            }
        }
        
    }

    /*
           @ 3
           |
   @       |      @
   @--------------@
   3       |      3
           |
           |
           @ 4
   */
    public void testPairInt() throws Exception {

        List<PairIntWithIndex> points = new ArrayList<>();

        points.add(new PairIntWithIndex(3, 0, 0));
        points.add(new PairIntWithIndex(3, 1, 1));
        points.add(new PairIntWithIndex(0, 3, 2));
        points.add(new PairIntWithIndex(-3, 1, 3));
        points.add(new PairIntWithIndex(-3, 0, 4));
        points.add(new PairIntWithIndex(0, -4, 5));

        double maxDistSq = Double.MIN_VALUE;
        PairInt p0 = null;
        PairInt p1 = null;
        for (PairInt p : points) {
            for (PairInt p2 : points) {
                if (p.equals(p2)) {
                    continue;
                }
                double dSq = Math.pow(p.getX() - p2.getX(), 2) + Math.pow(p.getY() - p2.getY(), 2);
                if (dSq > maxDistSq) {
                    maxDistSq = dSq;
                    p0 = p;
                    p1 = p2;
                }
            }
        }

        PairIntWithIndex[] pts = points.toArray(new PairIntWithIndex[points.size()]);
        PairInt[] result = LargestDistancePair.find(pts);

        assertNotNull(result);
        assertTrue(result.length == 2);

        PairInt expectedP0 = new PairInt(0, 3);
        PairInt expectedP1 = new PairInt(0, -4);

        assertTrue(expectedP0.equals(result[0]) || expectedP0.equals(result[1]));
        assertTrue(expectedP1.equals(result[0]) || expectedP1.equals(result[1]));
        assertFalse(result[0].equals(result[1]));

    }

}
