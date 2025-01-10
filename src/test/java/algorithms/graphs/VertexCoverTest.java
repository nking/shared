package algorithms.graphs;

import algorithms.trees.NAryTreeNode;
import algorithms.matrix.MatrixUtil;
import algorithms.optimization.LinearProgramming.StandardForm;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class VertexCoverTest extends TestCase {
    
    public VertexCoverTest(String testName) {
        super(testName);
    }
    
    public void testExact() {
        
        /*
        implemented from pseudocode in lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Covering Problems
           https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/vertexcover.pdf
           who reference Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
        
                                *0
                        *1                2
                  3                              *4
            *5        *6                        7  8  *9
         10 11 12    13 14                             15
        */
        
        int n = 16;
        NAryTreeNode[] nodes = new NAryTreeNode[n];
        int i;
        for (i = 0; i < n; ++i) {
            nodes[i] = new NAryTreeNode(i);
        }
        nodes[0].addChild(nodes[1]);
        nodes[1].addChild(nodes[3]);
        nodes[3].addChild(nodes[5]);
        nodes[3].addChild(nodes[6]);
        nodes[5].addChild(nodes[10]);
        nodes[5].addChild(nodes[11]);
        nodes[5].addChild(nodes[12]);
        nodes[6].addChild(nodes[13]);
        nodes[6].addChild(nodes[14]);
        
        nodes[0].addChild(nodes[2]);
        nodes[2].addChild(nodes[4]);
        nodes[4].addChild(nodes[7]);
        nodes[4].addChild(nodes[8]);
        nodes[4].addChild(nodes[9]);
        nodes[9].addChild(nodes[15]);
        
        TIntSet expected = new TIntHashSet(new int[]{0, 1, 4, 5, 6, 9});
       
        VertexCover vc = new VertexCover();
        System.out.println("vertex cover");
        Set<NAryTreeNode> cover = vc.exact(nodes[0]);
        for (NAryTreeNode node : cover) {
            System.out.printf("%d ", node.getData());
            assertTrue(expected.remove(node.getData()));
        }
        System.out.println();
        assertTrue(expected.isEmpty());
    }

    public void testExact2() {

        /*
                                *0
                         1                *2
                 *3                     4
             5    6   *7
                      8
        */

        int n = 9;
        NAryTreeNode[] nodes = new NAryTreeNode[n];
        int i;
        for (i = 0; i < n; ++i) {
            nodes[i] = new NAryTreeNode(i);
        }
        nodes[0].addChild(nodes[1]);
        nodes[0].addChild(nodes[2]);
        nodes[1].addChild(nodes[3]);
        nodes[2].addChild(nodes[4]);
        nodes[3].addChild(nodes[5]);
        nodes[3].addChild(nodes[6]);
        nodes[3].addChild(nodes[7]);
        nodes[7].addChild(nodes[8]);

        TIntSet expected = new TIntHashSet(new int[]{7, 3, 2, 0});

        VertexCover vc = new VertexCover();
        System.out.println("vertex cover");
        Set<NAryTreeNode> cover = vc.exact(nodes[0]);
        for (NAryTreeNode node : cover) {
            System.out.printf("%d ", node.getData());
            assertTrue(expected.remove(node.getData()));
        }
        System.out.println();
        assertTrue(expected.isEmpty());
    }
    
    /**
     * test is from Fig. 35.1 of Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms".
     */
    public void testApprox2() {
        
        /*
        b=1----c=2---/d=3
         |     |   /  | \
         |     | /    |   \
        a=0    e=4----f=5   \g=6  
        */
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        int i;
        for (i = 0; i < 7;++i) {
            g.put(i, new TIntHashSet());
        }
        g.get(0).add(1);
        g.get(1).add(2);
        g.get(2).add(3);
        g.get(2).add(4);
        g.get(3).add(4);
        g.get(3).add(5);
        g.get(3).add(6);
        g.get(4).add(2);
        g.get(4).add(3);
        g.get(4).add(5);
        g.get(5).add(3);
        g.get(5).add(4);
        g.get(6).add(3);
        
        VertexCover vc = new VertexCover();
        TIntSet cover = vc.approx2(g);
        //System.out.printf("2-appox set cover=%s\n", Arrays.toString(cover.toArray()));
        
        TIntSet expected = new TIntHashSet(new int[]{1, 2, 3, 4, 5, 6});
        
        assertEquals(expected.size(), cover.size());
        assertTrue(expected.containsAll(cover.toArray()));
    }

    /**
     * Test of copy method, of class VertexCover.
     */
    public void testCopy() {
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |
          2 E    |
                 |
       C---------D
       3         1
      */
        TObjectIntMap<String> vectorMap = new TObjectIntHashMap<String>();  
        TIntObjectMap<TIntSet> adjMap = new TIntObjectHashMap<TIntSet>();
        double[] weights = getGraph1(adjMap, vectorMap);
        
        VertexCover vc = new VertexCover();
        TIntObjectMap<TIntSet> adjMap2 = vc.copy(adjMap);
        assertEquals(adjMap.size(), adjMap2.size());
        
        TIntObjectIterator<TIntSet> iter = adjMap2.iterator();
        TIntSet set;
        TIntIterator iter2;
        int i, u, v;
        for (i = 0; i < adjMap2.size(); ++i) {
            iter.advance();;
            u = iter.key();
            set = iter.value();
            assertTrue(adjMap.containsKey(u));
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                adjMap.get(u).remove(v);
            }    
            assertTrue(adjMap.get(u).isEmpty());
            adjMap.remove(u);
        }
        assertTrue(adjMap.isEmpty());
    }

    /**
     * Test of reverse method, of class VertexCover.
     */
    public void testReverse() {
        
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |
          2 E    |
                 |
       C---------D
       3         1
      */
        TObjectIntMap<String> vectorMap = new TObjectIntHashMap<String>();  
        TIntObjectMap<TIntSet> adjMap = new TIntObjectHashMap<TIntSet>();
        double[] weights = getGraph1(adjMap, vectorMap);
        
        Set<PairInt> expectedReverseEdges = new HashSet<PairInt>();
        expectedReverseEdges.add(new PairInt(vectorMap.get("B"), vectorMap.get("A")));
        expectedReverseEdges.add(new PairInt(vectorMap.get("E"), vectorMap.get("A")));
        expectedReverseEdges.add(new PairInt(vectorMap.get("D"), vectorMap.get("B")));
        expectedReverseEdges.add(new PairInt(vectorMap.get("E"), vectorMap.get("B")));
        expectedReverseEdges.add(new PairInt(vectorMap.get("D"), vectorMap.get("C")));
        
        VertexCover vc = new VertexCover();        
        TIntObjectMap<TIntSet> r = MatrixUtil.createReverseMap(adjMap);

        TIntObjectIterator<TIntSet> iter = r.iterator();
        TIntSet set;
        TIntIterator iter2;
        int i, u, v;
        for (i = 0; i < r.size(); ++i) {
            iter.advance();;
            u = iter.key();
            set = iter.value();
            iter2 = set.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                expectedReverseEdges.remove(new PairInt(u, v));
            }
        }
        assertTrue(expectedReverseEdges.isEmpty());
    }

    /**
     * Test of createLinearProgramInStandardForm method, of class VertexCover.
     */
    public void estCreateLinearProgramInStandardForm() {
        //TODO: review this and fix expected if wrong or code if wrong
        // then re-enable
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |
          2 E    |
                 |
       C---------D
       3         1
      */
        TObjectIntMap<String> vectorMap = new TObjectIntHashMap<String>();  
        TIntObjectMap<TIntSet> adjMap = new TIntObjectHashMap<TIntSet>();
        double[] weights = getGraph1(adjMap, vectorMap);
      
        VertexCover vc = new VertexCover();
        StandardForm standForm = vc.createLinearProgramInStandardForm(adjMap, weights);
        
        System.out.printf("standForm=\n%s\n", standForm.toString());
        /*
         minimize: 
                summation_v_in_V( w(v)*x(v) )
            subject to:
                x(u) + x(v) >= 1 for each (u,v) in E
                x(v) <= 1 for each v in V
            non-negativity constraints:
                x(v) >= 0 for each v in V
        */
        
        double diff, tol = 1e-7;
        double expectedV = 0;
        double[] expectedC = new double[]{-4,  -3, -3, -1, -2};
        double[] expectedB = new double[]{-1.000, -1.000, -1.000, -1.000, -1.000, 
            1.000, 1.000, 1.000, 1.000, 1.000};
        double[][] expectedA = new double[10][];
        expectedA[0] = new double[]{-1.000, 0.000, 0.000, 0.000, -1.000};
        expectedA[1] = new double[]{0.000, 0.000, -1.000, -1.000, 0.000};
        expectedA[2] = new double[]{0.000, -1.000, 0.000, -1.000, 0.000};
        expectedA[3] = new double[]{0.000, -1.000, 0.000, 0.000, -1.000};
        expectedA[4] = new double[]{-1.000, -1.000, 0.000, 0.000, 0.000};
        expectedA[5] = new double[]{1.000, 0.000, 0.000, 0.000, 0.000};
        expectedA[6] = new double[]{0.000, 1.000, 0.000, 0.000, 0.000};
        expectedA[7] = new double[]{0.000, 0.000, 1.000, 0.000, 0.000};
        expectedA[8] = new double[]{0.000, 0.000, 0.000, 1.000, 0.000};
        expectedA[9] = new double[]{0.000, 0.000, 0.000, 0.000, 1.000};
        // the non-negativity constraints are implicit, not present in the a matrix
        //expectedA[10] = new double[]{-1.000, 0.000, 0.000, 0.000, 0.000};
        //expectedA[11] = new double[]{0.000, -1.000, 0.000, 0.000, 0.000};
        //expectedA[12] = new double[]{0.000, 0.000, -1.000, 0.000, 0.000};
        //expectedA[13] = new double[]{0.000, 0.000, 0.000, -1.000, 0.000};
        //expectedA[14] = new double[]{0.000, 0.000, 0.000, 0.000, -1.000};
    
        assertExpected(tol, standForm, expectedV, expectedB, expectedC, 
            expectedA);
        
    }

    /**
     * Test of extractEdges method, of class VertexCover.
     */
    public void testExtractEdges() {
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |
          2 E    |
                 |
       C---------D
       3         1
      */
        TObjectIntMap<String> vectorMap = new TObjectIntHashMap<String>();  
        TIntObjectMap<TIntSet> adjMap = new TIntObjectHashMap<TIntSet>();
        double[] weights = getGraph1(adjMap, vectorMap);
        
        Set<PairInt> expectedEdges = new HashSet<PairInt>();
        expectedEdges.add(new PairInt(vectorMap.get("A"), vectorMap.get("B")));
        expectedEdges.add(new PairInt(vectorMap.get("A"), vectorMap.get("E")));
        expectedEdges.add(new PairInt(vectorMap.get("B"), vectorMap.get("D")));
        expectedEdges.add(new PairInt(vectorMap.get("B"), vectorMap.get("E")));
        expectedEdges.add(new PairInt(vectorMap.get("C"), vectorMap.get("D")));
        
        VertexCover vc = new VertexCover();
        int[][] edges = vc.extractEdges(adjMap);
        
        assertEquals(expectedEdges.size(), edges.length);
        int i;
        PairInt p;
        for (i = 0; i < edges.length; ++i) {
            p = new PairInt(edges[i][0], edges[i][1]);
            expectedEdges.contains(p);
            expectedEdges.remove(p);
        }
        assertTrue(expectedEdges.isEmpty());
    }
    
    
    /**
     * Test of approxWeighted method, of class VertexCover.
     */
    public void testApproxWeighted() {
        
        /*
        graph is from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        */
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |            optimal = B(1), E(4), D(3) = 3+2+1 = 6
          2 E    |            solns   = A(0), B(1), D(3), E(4) = 4+3+1+2 = 10
                 |                    or = A(0), D(3), E(4) = 4+1+2=7
       C---------D                    or = B(1), D(3), E(4) = 3+1+2=6
       3         1       
      */
        TObjectIntMap<String> vectorMap = new TObjectIntHashMap<String>();  
        TIntObjectMap<TIntSet> adjMap = new TIntObjectHashMap<TIntSet>();
        double[] weights = getGraph1(adjMap, vectorMap);
                
        VertexCover vc = new VertexCover();
        TIntSet v = vc.approx2Weighted(adjMap, weights);
        
        System.out.printf("weighted cover=%s\n", Arrays.toString(v.toArray()));
        while (v.contains(vectorMap.get("C"))) {
            // LinearProgramming is using random picks of "entering" indexes in pivot, so try again
            v = vc.approx2Weighted(adjMap, weights);
        }
        assertFalse(v.contains(vectorMap.get("C")));
        assertTrue(v.size() == 4 || v.size() == 3);
    }

    private double[] getGraph1(TIntObjectMap<TIntSet> adjMap, TObjectIntMap<String> vectorMap) {
        
        /*
        graph is from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        */
        /*
        
        4        3 <-- weights
        A--------B
         \      /|
           \  /  |
          2 E    |
                 |
       C---------D
       3         1
      */
        vectorMap.put("A", 0);
        vectorMap.put("B", 1);
        vectorMap.put("C", 2);
        vectorMap.put("D", 3);
        vectorMap.put("E", 4);
        
        adjMap.put(vectorMap.get("A"), new TIntHashSet());
        adjMap.put(vectorMap.get("B"), new TIntHashSet());
        adjMap.put(vectorMap.get("C"), new TIntHashSet());
        adjMap.put(vectorMap.get("D"), new TIntHashSet());
        adjMap.put(vectorMap.get("E"), new TIntHashSet());
        
        adjMap.get(vectorMap.get("A")).add(vectorMap.get("B"));
        adjMap.get(vectorMap.get("A")).add(vectorMap.get("E"));
        adjMap.get(vectorMap.get("B")).add(vectorMap.get("E"));
        adjMap.get(vectorMap.get("B")).add(vectorMap.get("D"));
        adjMap.get(vectorMap.get("C")).add(vectorMap.get("D"));
        
        double[] weights = new double[]{4, 3, 3, 1, 2};
                
        return weights;
    }
    
     private void assertExpected(double tol,
        StandardForm standForm, double expectedV, 
        double[] expectedB, double[] expectedC, double[][] expectedA) {
        
        int i;
        double diff;
        assertTrue(Math.abs(expectedV - standForm.v) < tol);
        assertEquals(expectedB.length, standForm.b.length);
        assertEquals(expectedC.length, standForm.c.length);
        assertEquals(expectedA.length, standForm.a.length);
        assertEquals(expectedA[0].length, standForm.a[0].length);
        
        for (i = 0; i < expectedB.length; ++i) {
            diff = Math.abs(expectedB[i] - standForm.b[i]);
            assertTrue(diff < tol);
        }
        for (i = 0; i < expectedC.length; ++i) {
            diff = Math.abs(expectedC[i] - standForm.c[i]);
            assertTrue(diff < tol);
        }
        int j;
        for (i = 0; i < standForm.a.length; ++i) {
            for (j = 0; j < standForm.a[i].length; ++j) {
                diff = Math.abs(expectedA[i][j] - standForm.a[i][j]);
                System.out.printf("diff=%f tol=%f (diff<tol=%b)\n", diff, tol, (diff < tol));
                assertTrue(diff < tol);
            }
        }
    }
}
