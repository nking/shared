package algorithms.graphs;

import algorithms.graphs.DegeneracyOrderingMatulaBeck.TBucketQueue;
import algorithms.graphs.DegeneracyOrderingMatulaBeck.IBucketQueue;
import algorithms.graphs.DegeneracyOrderingMatulaBeck.SBucketQueue;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.LinkedHashMap;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DegeneracyOrderingMatulaBeckTest extends TestCase {
    
    public DegeneracyOrderingMatulaBeckTest(String testName) {
        super(testName);
    }
    
    public TIntObjectMap<TIntSet> getG0() {
        // from https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        int n = 6, i, j;
        for (i = 0; i < n; ++i) {
            g.put(i, new TIntHashSet());
        }
        g.get(1-1).add(5-1);
        g.get(1-1).add(2-1);
        g.get(2-1).add(5-1);
        g.get(2-1).add(3-1);
        g.get(3-1).add(4-1);
        g.get(3-1).add(2-1);
        g.get(4-1).add(6-1);
        g.get(4-1).add(5-1);
        g.get(5-1).add(4-1);
        g.get(5-1).add(2-1);
        g.get(5-1).add(1-1);
        g.get(6-1).add(4-1);
        return g;
    }
    
    public void testBucketQueue() {
        // from https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        
        TIntObjectMap<TIntSet> g = getG0();
        
        IBucketQueue bQ;
        
        for (int ii = 0; ii < 2; ++ii) {

            // -- test init and remove --
            if (ii == 0) {
                bQ = new TBucketQueue(g);
            } else {
                bQ = new SBucketQueue(g);
            }
            assertEquals(3, bQ.size());
            assertTrue(bQ.contains(1));
            assertTrue(bQ.contains(2));
            assertTrue(bQ.contains(3));
            assertFalse(bQ.isEmpty());

            int i, v;
            TIntSet bucket;
            while (bQ.size() > 0) {
                i = bQ.extractMinimumPt1();
                //System.out.printf("extractMin=%d: ", i);
                bucket = bQ.getBucket(i);
                while (!bucket.isEmpty()) {
                    v = bucket.iterator().next();
                    assertTrue(bQ.containsVertex(v));
                    //System.out.printf(" %d", v);
                    bQ.remove(i, v, bucket);
                }
                //System.out.printf("\n");

                // bucket for key i is empty, so key i should have been removed:
                assertFalse(bQ.contains(i));
            }
            //System.out.printf("bQ.size=%d\n", bQ.size());

            /*
            [junit] bM key=1, val=5, max=1
            [junit] bM key=3, val=4, max=3
            [junit] bM key=2, val=3, max=3
            [junit] bM key=2, val=2, max=3
            [junit] bM key=2, val=1, max=3
            [junit] bM key=2, val=0, max=3
            [junit] nBins=3  rt of ops=1.5849625007211563
            [junit] bQ key=3
            [junit] bQ key=2
            [junit] bQ key=1
            */

            // -- test moveItem --
            int ex = 0;
            try {
                bQ.moveItem(5, 1, 0);
            } catch(IllegalStateException e) {
                ex = 1;
            }
            assertEquals(1, ex);

            if (ii == 0) {
                bQ = new TBucketQueue(g);
            } else {
                bQ = new SBucketQueue(g);
            }
            bQ.moveItem(5, 1, 0);
            assertEquals(3, bQ.size());
            assertTrue(bQ.contains(0));
            assertFalse(bQ.contains(1));
            assertTrue(bQ.contains(2));
            assertTrue(bQ.contains(3));

            bQ.moveItem(5, 0, 2);
            assertEquals(2, bQ.size());
            assertTrue(bQ.contains(2));
            assertTrue(bQ.contains(3));
            assertFalse(bQ.contains(0));

            bQ.moveItem(5, 1);
            assertEquals(3, bQ.size());
            assertTrue(bQ.contains(2));
            assertTrue(bQ.contains(3));
            assertTrue(bQ.contains(1));
            assertFalse(bQ.contains(0));
            //
        }
    }

    public void testFindDegeneracyOrder() {
        // from https://en.m.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        
        TIntObjectMap<TIntSet> g = getG0();
        int n = g.size();
        
        int[] out = new int[n];
        
        int expResult = 0;
        
        int result = DegeneracyOrderingMatulaBeck.findDegeneracyOrder(g, out);
        
        System.out.printf("k=%d order=%s\n", result, Arrays.toString(out));
        /*
        maxDegree=3, n=6
        k=2 order=[2, 0, 1, 4, 3, 5]
                   3  1  2  5  4  6 <-- indexes in reference frame of wikipedia example

                   index 6 has degree 1 which is smallest, so remove it.
                   the graph at this point is k=2.
                   removing index 4 (d=2), then index 5 (d=2), index 2 (d=2), index 1 (d=1), index 3 (d=0)

         */
    }
    
}
