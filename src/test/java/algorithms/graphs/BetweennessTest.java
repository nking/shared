package algorithms.graphs;

import algorithms.graphs.Betweenness.Results;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BetweennessTest extends TestCase {
    
    public BetweennessTest(String testName) {
        super(testName);
    }
    
    public void testGirvanNewmanDAG() {
        int nV = 7;
        SimpleLinkedListNode[] adjList = new SimpleLinkedListNode[nV];
        for (int i = 0; i < nV; ++i) {
            adjList[i] = new SimpleLinkedListNode();
        }
        adjList[0].insert(1); adjList[0].insert(2);
        adjList[1].insert(3);
        adjList[2].insert(3); adjList[2].insert(4);
        adjList[3].insert(5);
        adjList[4].insert(5); adjList[4].insert(6);
        
        Betweenness b = new Betweenness();
        Results r;        
        TObjectFloatMap<PairInt> edges;
        TIntSet v;
        
        /* this shows that it doesn't matter which source index you choose,
        but would like to let the user to select it if they know something about the graph.
        */
        /*System.out.println("iterating over each as a source node:");
        for (int ii = 0; ii < nV; ++ii) {
            r = b.girvanNewmanDAG(adjList, ii);
            System.out.println("result for " + ii + " =\n" + r.toString());
        }*/
        
        TObjectFloatMap<PairInt> expected = new TObjectFloatHashMap<PairInt>();
        expected.put(new PairInt(0, 1), 11.f/6.f);
        expected.put(new PairInt(0, 2), 25.f/6.f);
        expected.put(new PairInt(1, 3), 5.f/6.f);
        expected.put(new PairInt(2, 3), 5.f/6.f);
        expected.put(new PairInt(3, 5), 2.f/3.f);
        expected.put(new PairInt(2, 4), 7.f/3.f);
        expected.put(new PairInt(4, 5), 1.f/3.f);
        expected.put(new PairInt(4, 6), 1.f);
        
        r = b.girvanNewmanDAG(adjList, 0);
        System.out.println("result=\n" + r.toString());
        
        edges = r.getEdges();
        assertNotNull(edges);
        assertEquals(expected.size(), edges.size());
        v = r.getVertexes();
        assertNotNull(v);
        assertEquals(nV, v.size());
        
        double tol = 1.e-3;
        TObjectFloatIterator<PairInt> iter = expected.iterator();
        
        for (int i = 0; i < nV; ++i) {
            assertTrue(v.contains(i));
            iter.advance();
            PairInt euv = iter.key();
            float ew = iter.value();
            assertTrue(edges.containsKey(euv));
            assertTrue(Math.abs(edges.get(euv) - ew) < tol);
        }
                
        // ------
        System.out.println("\ncutting edge=(0,2)");
        adjList[0].delete(2);
        
        r = b.girvanNewmanDAG(adjList, 0);
        edges = r.getEdges();
        
        System.out.println("result2=\n" + r.toString());
        expected = new TObjectFloatHashMap<PairInt>();
        expected.put(new PairInt(0, 1), 3.f);
        expected.put(new PairInt(1, 3), 2.f);
        expected.put(new PairInt(2, 3), 5.f/6.f);
        expected.put(new PairInt(3, 5), 2.f/3.f);
        expected.put(new PairInt(2, 4), 7.f/3.f);
        expected.put(new PairInt(4, 5), 1.f/3.f);
        expected.put(new PairInt(4, 6), 1.f);
        
        assertEquals(expected.size(), edges.size());
        assertEquals(nV, r.getVertexes().size());
        
        iter = expected.iterator();
        
        for (int i = 0; i < nV; ++i) {
            assertTrue(v.contains(i));
            iter.advance();
            PairInt euv = iter.key();
            float ew = iter.value();
            assertTrue(edges.containsKey(euv));
            assertTrue(Math.abs(edges.get(euv) - ew) < tol);
        }
    }
    
    public void testGirvanNewmanDirectionless() {
        int nV = 7;
        
        SimpleLinkedListNode[] adjList2 = new SimpleLinkedListNode[nV];
        for (int i = 0; i < nV; ++i) {
            adjList2[i] = new SimpleLinkedListNode();
        }
        adjList2[0].insert(1); adjList2[0].insert(2);
        adjList2[1].insert(3);
        adjList2[2].insert(3); adjList2[2].insert(4);
        adjList2[3].insert(5);
        adjList2[4].insert(5); adjList2[4].insert(6);
        
        adjList2[1].insert(0); adjList2[2].insert(0);
        adjList2[3].insert(1);
        adjList2[3].insert(2); adjList2[4].insert(2);
        adjList2[5].insert(3);
        adjList2[5].insert(4); adjList2[6].insert(4);
        
        Betweenness b = new Betweenness();
        Results r;        
        TObjectFloatMap<PairInt> edges;
        TIntSet v;
        
        TObjectFloatMap<PairInt> expected = new TObjectFloatHashMap<PairInt>();
        expected.put(new PairInt(0, 1), 11.f/6.f);
        expected.put(new PairInt(0, 2), 25.f/6.f);
        expected.put(new PairInt(1, 3), 5.f/6.f);
        expected.put(new PairInt(2, 3), 5.f/6.f);
        expected.put(new PairInt(3, 5), 2.f/3.f);
        expected.put(new PairInt(2, 4), 7.f/3.f);
        expected.put(new PairInt(4, 5), 1.f/3.f);
        expected.put(new PairInt(4, 6), 1.f);
        
        r = b.girvanNewmanDirectionless(adjList2, 0);
        System.out.println("result=\n" + r.toString());
        
        edges = r.getEdges();
        assertNotNull(edges);
        assertEquals(expected.size(), edges.size());
        v = r.getVertexes();
        assertNotNull(v);
        assertEquals(nV, v.size());
        
        double tol = 1.e-3;
        TObjectFloatIterator<PairInt> iter = expected.iterator();
        
        for (int i = 0; i < nV; ++i) {
            assertTrue(v.contains(i));
            iter.advance();
            PairInt euv = iter.key();
            float ew = iter.value();
            assertTrue(edges.containsKey(euv));
            assertTrue(Math.abs(edges.get(euv) - ew) < tol);
        }
        
    }
}
