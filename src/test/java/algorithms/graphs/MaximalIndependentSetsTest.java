package algorithms.graphs;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

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
