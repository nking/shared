package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PushbackInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * rough reduced function parser built to read the zachary's karate club file
 * from http://www-personal.umich.edu/~mejn/netdata/
 * into a graph.
 * 
 * TODO: should change the reading of characters to handle low and high surrogate code points
 * upon need.
 * 
 * @author nichole
 */
public class NewmanGMLParser {

    public static class GMLGraph {
        String graphType;
        TObjectFloatMap<PairInt> edgeWeightMap;
        TIntObjectMap<String> nodeIdLabelMap;
    }
        
    public static GMLGraph readGraph(String filePath) 
        throws FileNotFoundException, IOException {
        
        // format at
        //http://graphics.stanford.edu/courses/cs368-00-spring/TA/manuals/LEDA/gml_graph.html#29665
                
        // node map with key = id, value = label.
        TIntObjectMap<String> nodeIdLabelMap = new TIntObjectHashMap<String>();
        
        // edge map, key = (source, target), value = weight where source and target are edge node ids.
        TObjectFloatMap<PairInt> edgeMap = new TObjectFloatHashMap<PairInt>();
        
        String graphType = null;
        
        FileReader reader = null;
        BufferedReader in = null;
        
        // successfully read from stream
        boolean r = true;
        
        StringBuilder type = new StringBuilder();
        StringBuilder content = new StringBuilder();
        String parse = null;
        int eIdx, nIdx;
        
        try {
            
            in = new BufferedReader(new FileReader(new File(filePath)));

            // find graph [ which may have an end of line marker and spaces in between
            r = readToEndOf(in, "graph");
            r = readToEndOf(in, "[");
            
            // for first object read only , parse type for a graph type
            r = readTypeAndContent(in, type, content);
            parse = type.toString().trim();
            eIdx = parse.indexOf("node");
            nIdx = parse.indexOf("edge");
            if (eIdx != 0 && nIdx != 0) {
                if (eIdx > 0) {
                    graphType = parse.substring(0, eIdx);
                    //parseAndStoreEdge(content, edgeMap);
                } else if (nIdx > 0) {
                    graphType = parse.substring(0, nIdx);
                    //parseAndStoreNode(content, nodeIdLabelMap);
                } else {
                    throw new IllegalStateException("did not find any objects in graph");
                }
            }
            
            while (r) {
                r = readTypeAndContent(in, type, content);
                
                /*if (type.equals("node")){
                    parseAndStoreNode(content, nodeIdLabelMap);
                } else if (type.equals("edge")) {
                    parseAndStoreEdge(content, edgeMap);
                } else {
                    break;
                }*/
            }
            
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        GMLGraph g = new GMLGraph();
        g.edgeWeightMap = edgeMap;
        g.nodeIdLabelMap = nodeIdLabelMap;
        g.graphType = graphType;
        return g;
    }
    
    static boolean readToEndOf(BufferedReader in, String srch) throws IOException {
        int ch;
        //NOTE: for UTF-8 may have code characters that are high and low surrogate paris of code point
        //   so would need to change this
        while ((ch = in.read()) > -1) {
            srchloop :
            if ((char)ch == srch.charAt(0)) {
                for (int i = 1; i < srch.length(); ++i) {
                    ch = in.read();
                    if ((char)ch != srch.charAt(i)) {
                        break srchloop;
                    }
                }
                return true;
            }
        }
        return false;
    }
    
    /**
     * read everything up to first [ as type, then everything up until last ] as content
     * @param in
     * @param type
     * @param content
     * @return
     * @throws IOException 
     */
    static boolean readTypeAndContent(BufferedReader in, StringBuilder type, 
        StringBuilder content) throws IOException {
        int ch;
        //NOTE: for UTF-8 may have code characters that are high and low surrogate paris of code point
        //   so would need to change this
        while ((ch = in.read()) > -1) {
            //System.out.println("c=" + (char)ch);
            if ((char)ch == '[') {
                break;
            }
            type.append((char)ch);
        }
        while ((ch = in.read()) > -1) {
            //System.out.println("c=" + (char)ch);
            if ((char)ch == ']') {
                return true;
            }
            content.append((char)ch);
        }
        return false;
    }
    
}
