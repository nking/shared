package algorithms.graphs;

import algorithms.graphs.UnweightedGraphCommunityFinder.BestDivision;
import algorithms.util.ResourceFinder;
import algorithms.util.SimpleLinkedListNode;
import java.io.IOException;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class UnweightedGraphCommunityFinderTest extends TestCase {
    
    public UnweightedGraphCommunityFinderTest(String testName) {
        super(testName);
    }

    /**
     * Test of girvanNewman2002 method, of class UnweightedGraphCommunityFinder.
     */
    public void testGirvanNewman2002_0() throws IOException {
        
        String path = ResourceFinder.findFileInTestResources("karate.gml");
        
        NewmanGMLParser.GMLGraph g = NewmanGMLParser.readGraph(path);
        
        //TODO: change to method whch reads from indexes that start at "1" then
        //   remaps to "0" based graph.
        SimpleLinkedListNode[] adjList = GraphUtil.createAdjacencyList2(g);
        
        int src = 0;
        
        BestDivision best = UnweightedGraphCommunityFinder.girvanNewman2002(
            adjList, src);
        
    }
    
}
