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
    
    public void testMaximalWithBipartite() {
        
        System.out.println("testMaximalWithBipartite");
        
        TIntObjectMap<TIntSet> adj = getTestGraphCube();
        
        int nV = adj.size();
        
        // === build a cost matrix with all 1's excepting the diagonal and existing edges
        //     which are set to infinity ======
        
        float[][] matrix = new float[nV][];
        int i;
        int j;
        for (i = 0; i < nV; ++i) {
            matrix[i] = new float[nV];
            Arrays.fill(matrix[i], 1);
            matrix[i][i] = Float.POSITIVE_INFINITY;
        }
        TIntObjectIterator<TIntSet> iter = adj.iterator();
        TIntSet set;
        TIntIterator iter2;
        int u;
        int v;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            set = iter.value();
            if (set == null || set.isEmpty()) {
                continue;
            }
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                matrix[u][v] = Float.POSITIVE_INFINITY;
            }
        }
        
        //TOOD: consider cases where bipartite doesn't produce MIS.
        //      e.g. complete graphs?  every item in matrix would be infinity, so
        //       one would need to either to put responsiblity onto the user of graph not
        //       being complete, or assert the condition at start of method, or let
        //       bipartite proceed and then check that each matched pair is
        //       an independent set.  the later is useful for all graphs and
        //       only adds a linear runtime complexity to the algorithm.)
                        
        HungarianAlgorithm ha = new HungarianAlgorithm();
        int[][] matched = ha.computeAssignments(matrix);
        
        // === assert that each matched is an independent set and add it to a set of independent sets called s ======

        Set<TIntSet> s = new HashSet<TIntSet>();
        
        for (i = 0; i < matched.length; ++i) {
            System.out.printf("matched: %s\n", Arrays.toString(matched[i]));
            if (!(adj.containsKey(matched[i][0]) && adj.get(matched[i][0]).contains(matched[i][1])) 
              && !(adj.containsKey(matched[i][1]) && adj.get(matched[i][1]).contains(matched[i][0]))){
                s.add(new TIntHashSet(matched[i]));
            }
        }
        
        if (s.isEmpty()) {
            // return null
        }
        
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
        
        the maximal independent sets are
        {0, 4, 5, 7} <-- maximum
        {1, 2, 6, 3} <-- maximum
        {0, 6}
        {2, 7}
        {1, 5}
        {4, 3}
        */
        
        // find compatible
        boolean notAdded;
        boolean hasEdge;
        TIntSet m;
        int nMIS;
        List<TIntSet> mis = new ArrayList<TIntSet>();
        TIntIterator iter3;
        for (TIntSet si : s) {
            notAdded = true;
            nMIS = mis.size();
            for (j = 0; j < nMIS; ++j) {
                m = mis.get(j);
                // check whether any element of m, along with si elements, is an edge in graph G represented by adjMap
                hasEdge = false;
                iter2 = m.iterator();
                while (iter2.hasNext()) {
                    v = iter2.next();
                    iter3 = si.iterator();
                    while (!hasEdge && iter3.hasNext()) {
                        u = iter3.next();
                        if ((adj.containsKey(u) && adj.get(u).contains(v))
                        || (adj.containsKey(v) && adj.get(v).contains(u))) {
                            hasEdge = true;
                            break;
                        }
                    }
                }
                if (!hasEdge) {
                    m.addAll(si);
                    notAdded = false;
                }
            }
            if (notAdded) {
                mis.add(new TIntHashSet(si));
            }
        }
        
        for (TIntSet mi : mis) {
            System.out.printf("maximum: %s\n", Arrays.toString(mi.toArray()));
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
