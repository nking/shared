package algorithms.graphs;

import algorithms.graphs.Betweenness.Results;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
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
    
    public void testGirvanNewman() {
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
        Results r = b.girvanNewmanDistances(adjList, 0);
        
        TObjectFloatMap<PairInt> edges = r.getEdges();
        assertNotNull(edges);
        assertEquals(8, edges.size());
        TIntSet v = r.getVertexes();
        assertNotNull(v);
        assertEquals(7, v.size());
        assertEquals(0, r.getSrc());
        
        TObjectFloatMap<PairInt> expected = new TObjectFloatHashMap();
        expected.put(new PairInt(0, 1), 11.f/6.f);
        expected.put(new PairInt(0, 2), 25.f/6.f);
        expected.put(new PairInt(1, 3), 5.f/6.f);
        expected.put(new PairInt(2, 3), 5.f/6.f);
        expected.put(new PairInt(3, 5), 2.f/3.f);
        expected.put(new PairInt(2, 4), 7.f/3.f);
        expected.put(new PairInt(4, 5), 1.f/3.f);
        expected.put(new PairInt(4, 6), 1.f);
        
        for (int i = 0; i < nV; ++i) {
            assertTrue(v.contains(i));
        }
    }
}