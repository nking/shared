package algorithms.connected;

import algorithms.util.FormatArray;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinderTest extends TestCase {
    
    public ConnectedValuesGroupFinderTest(String testName) {
        super(testName);
    }
    
    public void test() {
        
        int w = 5;
        int h = 5;
       
        PixelHelper ph = new PixelHelper();
       
        List<TLongSet> expected = new ArrayList<TLongSet>();
        TLongSet a = new TLongHashSet();
        for (int i = 0; i < h; ++i) {
            a.add((0*w) + i);
            a.add((4*w) + i);
            a.add((i*w) + 0);
            a.add((i*w) + 4);
        }
        expected.add(a);
        a = new TLongHashSet();//1
        a.add((1*w) + 1);
        a.add((1*w) + 2);
        a.add((2*w) + 1);
        expected.add(a);
        a = new TLongHashSet();//2
        a.add((2*w) + 2);
        a.add((3*w) + 1);
        expected.add(a);
        a = new TLongHashSet();//3
        a.add((1*w) + 3);
        a.add((2*w) + 3);
        a.add((3*w) + 2);
        expected.add(a);
        
        int[][] data = new int[h][];
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};

        ConnectedValuesGroupFinder finder = new ConnectedValuesGroupFinder();
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
        
        //--------------------- test 2
        expected = new ArrayList<TLongSet>();
        a = new TLongHashSet();
        for (int i = 0; i < h; ++i) {
            a.add((0*w) + i);
            a.add((4*w) + i);
            a.add((i*w) + 0);
            a.add((i*w) + 4);
        }
        expected.add(a);
        a = new TLongHashSet();//1
        a.add((1*w) + 1);
        a.add((1*w) + 2);
        a.add((2*w) + 1);
        expected.add(a);
        a = new TLongHashSet();//3
        a.add((1*w) + 3);
        a.add((2*w) + 3);
        expected.add(a);
        /*
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        */
        finder = new ConnectedValuesGroupFinder();
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
        a.add((1*w) + 1);
        a.add((1*w) + 2);
        a.add((2*w) + 1);
        expected.add(a);
        a = new TLongHashSet();//3
        a.add((1*w) + 3);
        a.add((2*w) + 3);
        expected.add(a);
        /*
        data[0] = new int[]{0, 0, 0, 0, 0};
        data[1] = new int[]{0, 1, 1, 3, 0};
        data[2] = new int[]{0, 1, 2, 3, 0};
        data[3] = new int[]{0, 2, 3, 1, 0};
        data[4] = new int[]{0, 0, 0, 0, 0};
        */
        finder = new ConnectedValuesGroupFinder();
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
}
