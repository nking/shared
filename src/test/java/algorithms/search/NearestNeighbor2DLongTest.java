package algorithms.search;

import algorithms.misc.Misc0;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import static junit.framework.Assert.assertEquals;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NearestNeighbor2DLongTest extends TestCase {
    
    public NearestNeighbor2DLongTest() {
    }
    
    public void testLarge() {
        
        // exercising the code for large range dense filling
        
        //int w = 7000;
        //int h = 5000;
        //int bSz = 50;
        int w = 7000;
        int h = 500;
        int bSz = 50;
        
        TLongSet pixIdxs = new TLongHashSet();
        
        //randomly draw multiples of bSz in x and y and draw 
        //    squares that have random gaps added
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //seed = 1499996461259L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int nDraws = 25;
        int count = 0;
        
        int nBX = w/bSz;
        int nBY = h/bSz;
        Set<PairInt> chosen = new HashSet<PairInt>();
        int xc, yc;
        PairInt p;
        while (count < nDraws) {
            do {
                xc = rand.nextInt(nBX);
                yc = rand.nextInt(nBY);
                p = new PairInt(xc, yc);
            } while (chosen.contains(p));
            chosen.add(p);
            int x = xc*bSz;
            int y = yc*bSz;
            //System.out.format("%d:%d  %d:%d%n", x, x+bSz, y, y+bSz);
            draw(pixIdxs, w, h, x, y, bSz, rand);
            count++;
        }
        
        PixelHelper ph = new PixelHelper();
        
        long n2;
        //n2 = (w*h)/10;
        n2 = w*h/100;
        
        for (long i = 0; i < n2; ++i) {
            xc = rand.nextInt(w);
            yc = rand.nextInt(h);
            pixIdxs.add(ph.toPixelIndex(xc, yc, w));
        }
        
        Set<PairInt> nearest = null;
        
        NearestNeighbor2DLong nn2d = new NearestNeighbor2DLong(pixIdxs, w, h);
        nn2d.doNotUseCache();
        
        System.out.println("n in NN2D=" + pixIdxs.size() 
            + " nQueries=" + n2);
        
        for (long i = 0; i < n2; ++i) {
            xc = rand.nextInt(w);
            yc = rand.nextInt(h);
            
            nearest = nn2d.findClosest(xc, yc);
            nearest = nn2d.findClosest(xc, yc, bSz);
            nearest = nn2d.findClosestNotEqual(xc, yc);
            nearest = nn2d.findClosestWithinTolerance(xc, yc, 10.5);
        }
        
    }
    
    public void test0() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*        
          0  1  2  3  4  5  6  7  8  9  10
        0 0     2     4     6     8     
        1
        2 22    24    26    28    30    
        3 
        4 44    46          50    52    
        5
        6 66    68    70    72    74    
        7
        8 88    90    92    94    96    
        9
        */
        
        int maxX = 11;
        int maxY = 9;
        int k = 5;

        points.remove(new PairInt(4, 4));
        
        NearestNeighbor2DLong nn2D = new
            NearestNeighbor2DLong(points, maxX, maxY);
        
        Set<PairInt> nearest;
        
        nearest = nn2D.findClosest(2, 2);
        assertNotNull(nearest);
        assertTrue(nearest.iterator().next()
            .equals(new PairInt(2, 2)));
        
        nearest = nn2D.findClosest(4, 4);
                
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        assertEquals(expected.size(), nearest.size());
        
        for (PairInt p2 : nearest) {
            //System.out.println("p2=" + p2);
            assertTrue(expected.remove(p2));
        }
        
        assertEquals(0, expected.size());
        
        // -----------------------------------
        nearest = nn2D.findClosest(4, 4, 2);
        
        assertEquals(4, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
        // -----------------------------------
        points.add(new PairInt(4, 4));
        
        nn2D = new NearestNeighbor2DLong(points, maxX, maxY);
        
        nearest = nn2D.findClosest(4, 4);
        
        assertEquals(1, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
    }
    
    public void test0_wo_cache() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*        
          0  1  2  3  4  5  6  7  8  9  10
        0 0     2     4     6     8     
        1
        2 22    24    26    28    30    
        3 
        4 44    46          50    52    
        5
        6 66    68    70    72    74    
        7
        8 88    90    92    94    96    
        9
        */
        
        int maxX = 11;
        int maxY = 9;
        int k = 5;

        points.remove(new PairInt(4, 4));
        
        NearestNeighbor2DLong nn2D = new
            NearestNeighbor2DLong(points, maxX, maxY);
        
        nn2D.doNotUseCache();
        
        Set<PairInt> nearest;
        
        nearest = nn2D.findClosest(2, 2);
        assertNotNull(nearest);
        assertTrue(nearest.iterator().next()
            .equals(new PairInt(2, 2)));
        
        nearest = nn2D.findClosest(4, 4);
        
        assertEquals(4, nearest.size());
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
        // -----------------------------------
        nearest = nn2D.findClosest(4, 4, 2);
        
        assertEquals(4, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
        // -----------------------------------
        nearest = nn2D.findClosest(4, 4, 1);
        
        // -----------------------------------
        points.add(new PairInt(4, 4));
        
        nn2D = new NearestNeighbor2DLong(points, maxX, maxY);
        
        nearest = nn2D.findClosest(4, 4);
        
        assertEquals(1, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
    }
    
    public void test1() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*        
          0  1  2  3  4  5  6  7  8  9  10
        0 0     2     4     6     8     
        1
        2 22    24    26    28    30    
        3 
        4 44    46          50    52    
        5
        6 66    68    70    72    74    
        7
        8 88    90    92    94    96    
        9
        */
        
        int maxX = 11;
        int maxY = 9;
        int k = 5;
        
        NearestNeighbor2DLong nn2D = new
            NearestNeighbor2DLong(points, maxX, maxY);
        
        Set<PairInt> nearest = nn2D.findClosestWithinTolerance(
            4, 4, 2);
        
        assertEquals(5, nearest.size());
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
        expected.add(new PairInt(4, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
    }
    
    public void test1_pixelIdxs() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*        
          0  1  2  3  4  5  6  7  8  9  10
        0 0     2     4     6     8     
        1
        2 22    24    26    28    30    
        3 
        4 44    46    48    50    52    
        5
        6 66    68    70    72    74    
        7
        8 88    90    92    94    96    
        9
        */
        
        int maxX = 11;
        int maxY = 9;
        int k = 5;
        
        PixelHelper ph = new PixelHelper();
        TLongSet pixelIdxs = new TLongHashSet();
        for (PairInt p : points) {
            long pixIdx = ph.toPixelIndex(p, maxX);
            pixelIdxs.add(pixIdx);
        }
        
        NearestNeighbor2DLong nn2D = new
            NearestNeighbor2DLong(pixelIdxs, maxX, maxY);
        
        Set<PairInt> nearest = nn2D.findClosestWithinTolerance(
            4, 4, 2);
        
        assertEquals(5, nearest.size());
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
        expected.add(new PairInt(4, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
    }
    
    public void test1_NE() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*        
          0  1  2  3  4  5  6  7  8  9  10
        0 0     2     4     6     8     
        1
        2 22    24    26    28    30    
        3 
        4 44    46          50    52    
        5
        6 66    68    70    72    74    
        7
        8 88    90    92    94    96    
        9
        */
        
        int maxX = 11;
        int maxY = 9;
        int k = 5;
        
        PixelHelper ph = new PixelHelper();
        TLongSet pixelIdxs = new TLongHashSet();
        for (PairInt p : points) {
            long pixIdx = ph.toPixelIndex(p, maxX);
            pixelIdxs.add(pixIdx);
        }
        
        NearestNeighbor2DLong nn2D = new
            NearestNeighbor2DLong(pixelIdxs, maxX, maxY);
        
        Set<PairInt> nearest = nn2D.findClosestNotEqual(
            4, 4);
                
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
        
        assertEquals(expected.size(), nearest.size());
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
    }
    
    private Set<PairInt> getTestData() {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < 10; i += 2) {
            for (int j = 0; j < 10; j += 2) {
                PairInt p = new PairInt(i, j);
                points.add(p);
            }
        }
        
        return points;
    }
 
    private void draw(TLongSet pixIdxs, int w, int h, 
        int xc, int yc, int bSz, Random rand) {

        PixelHelper ph = new PixelHelper();
        
        for (int i = xc; i < xc + bSz; ++i) {
            if (i >= w) {
                break;
            }
            for (int j = yc; j < yc + bSz; ++j) {
                if (j >= h) {
                    break;
                }
                pixIdxs.add(ph.toPixelIndex(i, j, w));
                j += rand.nextInt(bSz/2);
            }
            i += rand.nextInt(bSz/2);
        }
    }
}
