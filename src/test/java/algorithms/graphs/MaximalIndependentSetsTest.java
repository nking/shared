package algorithms.graphs;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class MaximalIndependentSetsTest extends TestCase {
    
    public MaximalIndependentSetsTest(String testName) {
        super(testName);
    }

    public void testFindOneProbabilistic() {
        System.out.println("findOneProbabilistic");
        
        TIntObjectMap<TIntSet> adjMap = getTestGraphCube();
        
        // result is expected to be 1 within this list:
        List<TIntSet> expResult = getTestGraphCubeExpectedMIS();
        
        TIntSet result = MaximalIndependentSets.findOneProbabilistic(adjMap);
        System.out.printf("result=%s\n", Arrays.toString(result.toArray()));
        
        boolean found = false;
        
        int[] ex;
        int[] r = result.toArray();
        Arrays.sort(r);
        for (int i = 0; i < expResult.size(); ++i) {
            ex = expResult.get(i).toArray();
            Arrays.sort(ex);
            if (Arrays.equals(ex, r)) {
                found = true;
                break;
            }
        }
        assertTrue(found);
    }

   
    public void testFindOne() {
        System.out.println("findOne");
        
        TIntObjectMap<TIntSet> adjMap = getTestGraphCube();
        
        // result is expected to be 1 within this list:
        List<TIntSet> expResult = getTestGraphCubeExpectedMIS();
        
        TIntSet result = MaximalIndependentSets.findOne(adjMap);
        System.out.printf("result=%s\n", Arrays.toString(result.toArray()));
        
        boolean found = false;
        
        int[] ex;
        int[] r = result.toArray();
        Arrays.sort(r);
        for (int i = 0; i < expResult.size(); ++i) {
            ex = expResult.get(i).toArray();
            Arrays.sort(ex);
            if (Arrays.equals(ex, r)) {
                found = true;
                break;
            }
        }
        assertTrue(found);
    }
    
    public void testMaximumWithBipartite() {
        
        System.out.println("testMaximumWithBipartite");
        
        TIntObjectMap<TIntSet> adj = getTestGraphCube();
        
        /*
        matched: [0, 4]
        matched: [1, 2]
        matched: [2, 1]
        matched: [3, 6]
        matched: [4, 0]
        matched: [5, 7]
        matched: [6, 3]
        matched: [7, 5]
        
        s:
        {0, 4}
        {1, 2}
        {3, 6}
        {5, 7}
        
        compatible sets:
        iterate over each set to add to existing or start new
        out = compatible sets (no 2 vertices in any set are an edge in g)
        
        i  s[i]    out
        0  {0, 4}  [{0,4}]
        1  {1, 2}  [{0,4}],[{1, 2}]
        2  {3, 6}  [{0,4}],[{1, 2}, {3, 6}]
        3  {5, 7}  [{0,4},{5, 7}],[{1, 2}, {3, 6}]
        
        runtime: 2*1 + 2*2 + 2*4 + 2*6 = 2 + 24 = 2 + summation(2*(2*i) from i=1 to n-1)
        so is O(n)
        
        the maximum maximal independent sets are
        {0, 4, 5, 7} <-- maximum
        {1, 2, 3, 6} <-- maximum
        {0, 6}
        {2, 7}
        {1, 5}
        {4, 3}
        */
        
        List<TIntSet> expectedM = new ArrayList<TIntSet>();
        expectedM.add(new TIntHashSet(new int[]{0, 4, 5, 7}));
        expectedM.add(new TIntHashSet(new int[]{1, 2, 3, 6}));
        
        List<TIntSet> mis = MaximalIndependentSets.findAllMaximum(adj, adj.size());
        
        assertEquals(expectedM.size(), mis.size());
        
        int[] o1;
        int[] o2;
        for (TIntSet mi : mis) {
            System.out.printf("maximum: %s\n", Arrays.toString(mi.toArray()));
            assertTrue(expectedM.contains(mi));
        }
    }

    public static TIntObjectMap<TIntSet> getTestGraphCube() {
        // from https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        
        TIntObjectMap<TIntSet> adj = new TIntObjectHashMap<TIntSet>();
        int n = 8;
        int i;
        for (i = 0; i < n; ++i) {
            adj.put(i, new TIntHashSet());
        }
        
        adj.get(0).add(1);
        adj.get(0).add(2);
        adj.get(0).add(3);
        adj.get(1).add(0);
        adj.get(1).add(4);
        adj.get(1).add(7);
        adj.get(2).add(0);
        adj.get(2).add(4);
        adj.get(2).add(5);
        adj.get(3).add(0);
        adj.get(3).add(5);
        adj.get(3).add(7);
        adj.get(4).add(1);
        adj.get(4).add(2);
        adj.get(4).add(6);
        adj.get(5).add(2);
        adj.get(5).add(3);
        adj.get(5).add(6);
        adj.get(6).add(4);
        adj.get(6).add(5);
        adj.get(6).add(7);
        adj.get(7).add(1);
        adj.get(7).add(3);
        adj.get(7).add(6);
        
        return adj;
    }
    
    
    public static List<TIntSet> getTestGraphCubeExpectedMIS() {
        // from https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        List<TIntSet> out = new ArrayList<TIntSet>();
        
        out.add(new TIntHashSet(new int[]{0, 6}));
        out.add(new TIntHashSet(new int[]{0, 4, 5, 7}));
        out.add(new TIntHashSet(new int[]{2, 7}));
        out.add(new TIntHashSet(new int[]{1, 5}));
        out.add(new TIntHashSet(new int[]{1, 2, 6, 3}));
        out.add(new TIntHashSet(new int[]{4, 3}));
        
        return out;
    }
    
    public static List<TIntSet> getTestGraphCubeExpectedMaximum() {
        // from https://en.wikipedia.org/wiki/Maximal_independent_set#Listing_all_maximal_independent_sets
        List<TIntSet> out = new ArrayList<TIntSet>();
        
        out.add(new TIntHashSet(new int[]{0, 4, 5, 7}));
        out.add(new TIntHashSet(new int[]{1, 2, 6, 3}));
        
        return out;
    }
}
