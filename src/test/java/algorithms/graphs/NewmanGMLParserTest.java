package algorithms.graphs;

import algorithms.graphs.NewmanGMLParser.GMLGraph;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.list.TFloatList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NewmanGMLParserTest extends TestCase {
    
    public NewmanGMLParserTest(String testName) {
        super(testName);
    }

    /**
     * Test of readGraph method, of class NewmanGMLParser.
     */
    public void testReadGraph0() throws Exception {
        
        String path = ResourceFinder.findFileInTestResources("celegansneural.gml");
        GMLGraph g = NewmanGMLParser.readGraph(path);
        
        assertNotNull(g.graphType);
        assertEquals("directed 1", g.graphType);
        
        assertEquals(297, g.nodeIdLabelMap.size());
        
        // 2345 unique (u, v) edge keys, and several have more than one edge
        //    between the same pair of vertices.
        int nEExpected = 2359;
        int nE = 0;
        Map<PairInt, TFloatList> edgeMap = g.edgeWeightMap;
        Iterator<Entry<PairInt, TFloatList>> iter = edgeMap.entrySet().iterator();
        Entry<PairInt, TFloatList> p;
        while (iter.hasNext()) {
            p = iter.next();
            assertNotNull(p.getValue());
            nE += p.getValue().size();
        }
        assertEquals(nEExpected, nE);
        
        SimpleLinkedListNode[] a = GraphUtil.createAdjacencyList(g);
        assertEquals(297, a.length);
        
        iter = edgeMap.entrySet().iterator();
        PairInt uv;
        SimpleLinkedListNode node;
        while (iter.hasNext()) {
            p = iter.next();
            uv = p.getKey();
            node = a[uv.getX()];
            assertTrue(node.contains(uv.getY()));
        }
        
    }
    
    public void testReadGraph1() throws Exception {
        
        // indexes start at "1"
        String path = ResourceFinder.findFileInTestResources("karate.gml");
        GMLGraph g = NewmanGMLParser.readGraph(path);
        
        assertEquals(34, g.nodeIdLabelMap.size());
        
        int nEExpected = 78;
        int nE = 0;
        Map<PairInt, TFloatList> edgeMap = g.edgeWeightMap;
        Iterator<Entry<PairInt, TFloatList>> iter = edgeMap.entrySet().iterator();
        Entry<PairInt, TFloatList> p;
        while (iter.hasNext()) {
            p = iter.next();
            assertNotNull(p.getValue());
            nE += p.getValue().size();
        }
        assertEquals(nEExpected, nE);
        
        SimpleLinkedListNode[] a = GraphUtil.createAdjacencyList2(g);
        assertEquals(34 + 1, a.length);
        
        iter = edgeMap.entrySet().iterator();
        PairInt uv;
        SimpleLinkedListNode node;
        while (iter.hasNext()) {
            p = iter.next();
            uv = p.getKey();
            node = a[uv.getX()];
            assertTrue(node.contains(uv.getY()));
        }
    }

    /**
     * Test of readToEndOf method, of class NewmanGMLParser.
     */
    public void testReadToEndOf() throws Exception {
        String path = ResourceFinder.findFileInTestResources("celegansneural.gml");
        
        BufferedReader in = new BufferedReader(new FileReader(new File(path)));
            
        boolean a = NewmanGMLParser.readToEndOf(in, "graph");
        assertTrue(a);
        a = NewmanGMLParser.readToEndOf(in, "[");
        assertTrue(a);
        
        StringBuilder type = new StringBuilder();
        StringBuilder content = new StringBuilder();
        
        a = NewmanGMLParser.readTypeAndContent(in, type, content);
        //System.out.println("type=" + type.toString());
        //System.out.println("content=" + content.toString());
        
        assertTrue(a);
        
        in.close();
    }
    
}
