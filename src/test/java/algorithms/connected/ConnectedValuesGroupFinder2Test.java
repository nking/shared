package algorithms.connected;

import algorithms.misc.Misc0;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinder2Test extends TestCase {
    
    public ConnectedValuesGroupFinder2Test(String testName) {
        super(testName);
    }
    
    public void testLarge() {
        
        int w = 7000;
        int h = 5000;
        int bSz = 50;
        
        int c = 0;
        int[][] data = new int[w][];
        for (int i = 0; i < w; ++i) {
            data[i] = new int[h];
        }
        
        //randomly draw multiples of bSz in x and y and draw 
        //    increasing value squares.
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //seed = 1499572125948L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        int nDraws = 25;
        int count = 0;
        
        int nBX = w/bSz;
        int nBY = h/bSz;
        int v = 10;
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
            draw(data, x, y, bSz, v);
            count++;
            v *= 2;
        }
        
        ConnectedValuesGroupFinder2 finder = new ConnectedValuesGroupFinder2();
        finder.setMinimumNumberInCluster(2);
        finder.setToUse8Neighbors();
        TIntSet exclude = new TIntHashSet();
        exclude.add(0);
        finder.setValuesToExclude(exclude);
        List<TLongSet> groups = finder.findGroups(data);
        assertEquals(nDraws, groups.size());        
    }
    
    public void test() {
        
        int w = 5;
        int h = 5;
       
        PixelHelper ph = new PixelHelper();
       
        List<TLongSet> expected = new ArrayList<TLongSet>();
        TLongSet a = new TLongHashSet();
        for (int i = 0; i < h; ++i) {
            a.add(ph.toPixelIndex(0, i, w));
            a.add(ph.toPixelIndex(4, i, w));
            a.add(ph.toPixelIndex(i, 0, w));
            a.add(ph.toPixelIndex(i, 4, w));
        }
        expected.add(a);
        assertEquals(16, a.size());
        a = new TLongHashSet();//1
        a.add(ph.toPixelIndex(1, 1, w));
        a.add(ph.toPixelIndex(1, 2, w));
        a.add(ph.toPixelIndex(2, 1, w));
        expected.add(a);
        a = new TLongHashSet();//2
        a.add(ph.toPixelIndex(2, 2, w));
        a.add(ph.toPixelIndex(3, 1, w));
        expected.add(a);
        a = new TLongHashSet();//3
        a.add(ph.toPixelIndex(1, 3, w));
        a.add(ph.toPixelIndex(2, 3, w));
        a.add(ph.toPixelIndex(3, 2, w));
        expected.add(a);
        
        int[][] data = new int[w][];
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        
        ConnectedValuesGroupFinder2 finder = new ConnectedValuesGroupFinder2();
        finder.setMinimumNumberInCluster(2);
        finder.setToUse8Neighbors();
        List<TLongSet> groups = finder.findGroups(data);
        assertEquals(4, groups.size());
        assertEquals(4, expected.size());
        
        for (int i = 0; i < groups.size(); ++i) {
            TLongSet group = groups.get(i);
            
            long t = group.iterator().next();
            TLongSet e = null;
            for (int j = 0; j < expected.size(); ++j) {
                TLongSet e2 = expected.get(j);
                if (e2.contains(t)) {
                    e = e2; break;
                }
            }
            assertNotNull(e);
            
            TLongIterator iter = group.iterator();
            while (iter.hasNext()) {
                long pIdx = iter.next();
                assertTrue(e.remove(pIdx));
            }
            assertTrue(e.isEmpty());
            expected.remove(e);
        }
        
        //-------
        expected = new ArrayList<TLongSet>();
        a = new TLongHashSet();
        for (int i = 0; i < h; ++i) {
            a.add(ph.toPixelIndex(0, i, w));
            a.add(ph.toPixelIndex(4, i, w));
            a.add(ph.toPixelIndex(i, 0, w));
            a.add(ph.toPixelIndex(i, 4, w));
        }
        expected.add(a);
        assertEquals(16, a.size());
        a = new TLongHashSet();//1
        a.add(ph.toPixelIndex(1, 1, w));
        a.add(ph.toPixelIndex(1, 2, w));
        a.add(ph.toPixelIndex(2, 1, w));
        expected.add(a);
        a = new TLongHashSet();//3
        a.add(ph.toPixelIndex(1, 3, w));
        a.add(ph.toPixelIndex(2, 3, w));
        expected.add(a);
        /*
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        */
        finder = new ConnectedValuesGroupFinder2();
        finder.setMinimumNumberInCluster(2);
        assertTrue(finder.use4Neighbors);
        
        groups = finder.findGroups(data);
        
        assertEquals(3, groups.size());
    
        assertEquals(expected.size(), groups.size());
        
        for (int i = 0; i < groups.size(); ++i) {
            TLongSet group = groups.get(i);
            
            long t = group.iterator().next();
            TLongSet e = null;
            for (int j = 0; j < expected.size(); ++j) {
                TLongSet e2 = expected.get(j);
                if (e2.contains(t)) {
                    e = e2; break;
                }
            }
            assertNotNull(e);
            
            TLongIterator iter = group.iterator();
            while (iter.hasNext()) {
                long pIdx = iter.next();
                assertTrue(e.remove(pIdx));
            }
            assertTrue(e.isEmpty());
            expected.remove(e);
        }
        
         //-------
        expected = new ArrayList<TLongSet>();
        a = new TLongHashSet();//1
        a.add(ph.toPixelIndex(1, 1, w));
        a.add(ph.toPixelIndex(1, 2, w));
        a.add(ph.toPixelIndex(2, 1, w));
        expected.add(a);
        a = new TLongHashSet();//3
        a.add(ph.toPixelIndex(1, 3, w));
        a.add(ph.toPixelIndex(2, 3, w));
        expected.add(a);
        /*
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        */
        finder = new ConnectedValuesGroupFinder2();
        TIntSet exclude = new TIntHashSet();
        exclude.add(0);
        finder.setValuesToExclude(exclude);
        finder.setMinimumNumberInCluster(2);
        assertTrue(finder.use4Neighbors);
        
        groups = finder.findGroups(data);
            
        assertEquals(expected.size(), groups.size());
        
        for (int i = 0; i < groups.size(); ++i) {
            TLongSet group = groups.get(i);
            
            long t = group.iterator().next();
            TLongSet e = null;
            for (int j = 0; j < expected.size(); ++j) {
                TLongSet e2 = expected.get(j);
                if (e2.contains(t)) {
                    e = e2; break;
                }
            }
            assertNotNull(e);
            
            TLongIterator iter = group.iterator();
            while (iter.hasNext()) {
                long pIdx = iter.next();
                assertTrue(e.remove(pIdx));
            }
            assertTrue(e.isEmpty());
            expected.remove(e);
        }
    }

    private void draw(int[][] data, int xc, int yc, int bSz, int v) {

        for (int i = xc; i < xc + bSz; ++i) {
            if (i >= data.length) {
                break;
            }
            for (int j = yc; j < yc + bSz; ++j) {
                if (j >= data[i].length) {
                    break;
                }
                data[i][j] = v;
            }
        }
    }
}
