package algorithms.graphs;

import algorithms.graphs.NewmanGMLParser.GMLGraph;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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
    public void testReadGraph() throws Exception {
        
        String path = ResourceFinder.findFileInTestResources("celegansneural.gml");
        GMLGraph g = NewmanGMLParser.readGraph(path);
        assertNotNull(g.graphType);
        System.out.println("graphType=" + g.graphType);
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
        System.out.println("type=" + type.toString());
        System.out.println("content=" + content.toString());
        
        assertTrue(a);
        
        in.close();
    }
    
}
