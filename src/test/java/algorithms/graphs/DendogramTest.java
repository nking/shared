package algorithms.graphs;

import algorithms.util.ResourceFinder;
import algorithms.util.SimpleLinkedListNode;
import java.io.IOException;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DendogramTest extends TestCase {
    
    public DendogramTest(String testName) {
        super(testName);
    }
    
    public void test0() throws IOException {
        
        String path = ResourceFinder.findFileInTestResources("karate.gml");
        
        NewmanGMLParser.GMLGraph g = NewmanGMLParser.readGraph(path);
        
        SimpleLinkedListNode[] adjList = GraphUtil.createAdjacencyList2(g);
        
        Dendogram d = new Dendogram(adjList);
        
        int kFinal = d.createUsingGirvanNewman(0);
        
        List<Dendogram.DendogramLayer> layers = d.getLayers();
    }
}
