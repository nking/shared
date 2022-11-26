package algorithms.graphs;

import algorithms.util.PairInt;
import gnu.trove.list.TFloatList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * rough reduced function parser built to read a few of the GML files
 * from http://www-personal.umich.edu/~mejn/netdata/
 * into a graph.
 * 
 * the GML format is from:
   "GML: A portable Graph File Format" by Michael Himsolt
    http://openmis.ru/doc/clang/gml-tr.html 
 
 NOTE: if one wanted to define all required and optional attributes, could
 make a more formal parser using ANTLR at antlr.org
 
 * @author nichole
 */
public class NewmanGMLParser {

    private static Logger log = Logger.getLogger("NewmanGMLParser");

    public static class GMLGraph {
        String graphType;
        Map<PairInt, TFloatList> edgeWeightMap;
        TIntObjectMap<String> nodeIdLabelMap;
    }
        
    public static GMLGraph readGraph(String filePath) 
        throws FileNotFoundException, IOException {
        
        // format at
        //http://graphics.stanford.edu/courses/cs368-00-spring/TA/manuals/LEDA/gml_graph.html#29665
                
        // node map with key = id, value = label.
        TIntObjectMap<String> nodeIdLabelMap = new TIntObjectHashMap<String>();
        
        // edge map, key = (source, target), value = weights;  source and target are edge node ids.
        Map<PairInt, TFloatList> edgeMap = new HashMap<PairInt, TFloatList>();
        
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
            reader = new FileReader(new File(filePath));
            in = new BufferedReader(reader);

            // find graph [ which may have an end of line marker and spaces in between
            r = readToEndOf(in, "graph");
            r = readToEndOf(in, "[");
            
            // for first object read only , parse type for a graph type
            r = readTypeAndContent(in, type, content);
            parse = type.toString().trim();
            nIdx = parse.indexOf("node");
            eIdx = parse.indexOf("edge");
            if (eIdx != -1) {
                graphType = parse.substring(0, eIdx).trim();
                parseAndStoreEdge(content.toString(), edgeMap);
            } else if (nIdx != -1) {
                graphType = parse.substring(0, nIdx).trim();
                parseAndStoreNode(content.toString(), nodeIdLabelMap);
            } else {
                //empty graph
                return new GMLGraph();
            }
            
            while (r) {
                content = content.delete(0, content.length());
                type = type.delete(0, type.length());
                
                r = readTypeAndContent(in, type, content);
                parse = type.toString().trim();
                //System.out.printf(" %s", parse);   
                
                if (parse.equals("node")){
                    parseAndStoreNode(content.toString(), nodeIdLabelMap);
                } else if (parse.equals("edge")) {
                    parseAndStoreEdge(content.toString(), edgeMap);
                } else {
                    break;
                }
            }
            
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (in != null) {
                in.close();
            }
        }

        GMLGraph g = new GMLGraph();
        g.edgeWeightMap = edgeMap;
        g.nodeIdLabelMap = nodeIdLabelMap;
        g.graphType = graphType;
        return g;
    }
    
    /**
     * note: assumes srch characters each fit into 1 code point.
     * @param in
     * @param srch
     * @return
     * @throws IOException 
     */
    static boolean readToEndOf(BufferedReader in, String srch) throws IOException {
        int ch;
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
        while ((ch = in.read()) > -1) {
            if ((char)ch == '[') {
                break;
            }
            type.appendCodePoint(ch);
        }
        while ((ch = in.read()) > -1) {
            if ((char)ch == ']') {
                return true;
            }
            content.appendCodePoint(ch);
        }
        return false;
    }
    
    /**
     * parses the key-value pairs of the edge content.
     * recognizes source, target, and value.
     * the value by default is "1" if not present.
     * NOTE: the other attributes such as label are not currently stored.
     * 
     * @param content
     * @param edgeMap 
     */
    private static void parseAndStoreEdge(String content, Map<PairInt, TFloatList> edgeMap) {
        
        content = content.trim();
        
        String[] p = content.split("\\s+");
        //System.out.println(Arrays.toString(p));
        
        //source, target, and value
        String source = null;
        String target = null;
        String value = null;
        
        for (int i = 0; i < p.length; i+=2) {
            if (p[i].equalsIgnoreCase("source")) {
                source = p[i+1];
            } else if (p[i].equalsIgnoreCase("target")) {
                target = p[i+1];
            } else if (p[i].equalsIgnoreCase("value")) {
                value = p[i+1];
            }
        }
        
        if (source == null || target == null) {
            throw new IllegalArgumentException("edge is missing source and/or target: " + content);
        }

        PairInt e = new PairInt(Integer.parseInt(source), Integer.parseInt(target));
        
        if (!edgeMap.containsKey(e)) {
            edgeMap.put(e, new TFloatArrayList());
        } else {
            log.warning("edge " + e.toString() + " is in file more than once. ");
        }
        
        if (value == null) {
            edgeMap.get(e).add(1);
        } else {
            edgeMap.get(e).add(Float.parseFloat(value));
        }
    }
    
    /**
     * parses the key-value pairs from the node content and stores the
     * attributes id and label.  if label is not present, a string
     * representation of id's value is stored for it.
     * 
     * @param content
     * @param nodeIdLabelMap 
     */
    private static void parseAndStoreNode(String content, 
        TIntObjectMap<String> nodeIdLabelMap) {
        
        content = content.trim();
        
        String[] p = content.split("\\s+");
        
        String id = null;
        String label = null;
        
        for (int i = 0; i < p.length; i+=2) {
            if (p[i].equalsIgnoreCase("id")) {
                id = p[i+1];
            } else if (p[i].equalsIgnoreCase("label")) {
                label = p[i+1];
                label = label.replaceAll("^\"", "");
                label = label.replaceAll("\"$", "");
            }
        } 
        
        if (label == null && id != null) {
            label = id;
        }
        
        if (id == null) {
            throw new IllegalArgumentException("node is missing id: " + content);
        }
        
        int idI = Integer.parseInt(id);
        
        if (nodeIdLabelMap.containsKey(idI)) {
            log.warning("node id=" + id + " is in file more than once");
        }
        
        nodeIdLabelMap.put(idI, label);
    }

    /*
    from: "GML: A portable Graph File Format" by Michael Himsolt
    http://openmis.ru/doc/clang/gml-tr.html
    
    Global Defined Keys
        .id int
        Defines an identification number for an object. This is usually used to represent pointers.
        .label string
        Defines a label attached to an object.
        .comment string
        Defines a comment embedded in a GML file. Comments are ignored by the application.
        .Creator string
        Shows which application created this file and should therefore only be used once per file at the top level. .Creator is obviously unsafe.
        .graphics list
        Describes graphics which are used to draw a particular object.Within graphics, the following keys are defined:
        .graphics.x float
        Defines the x coordinate of the center of the object.
        .graphics.y float
        Defines the y coordinate of the center of the object.
        .graphics.z float
        Defines the z coordinate of the center of the object.
        .graphics.w float
        Defines the width of the object.
        .graphics.h float
        Defines the height of the object.
        .graphics.d float
        Defines the depth of the object.
        Coordinates are pixel coordinates on a standard 72 dpi drawing area. Applications may use them as screen coordinates.

    Example:
        graph [
          node [
            id 7
            label "5"
            edgeAnchor "corners"
            labelAnchor "n"
            graphics [
              center [ x 82.0000 y 42.0000 ]
              w 16.0000
              h 16.0000
              type "rectangle"
              fill "#000000"
            ]
          ]
          node [
            id 15
            label "13"
            edgeAnchor "corners"
            labelAnchor "c"
            graphics [
              center [ x 73.0000 y 160.000 ]
              w 16.0000
              h 16.0000
              type "rectangle"
              fill "#FF0000"
            ]
          ]
          edge [
            label "24"
            labelAnchor "first"
            source 7
            target 15
            graphics [
              type "line"
              arrow "last"
              Line [
                point [ x 82.0000 y 42.0000 ]
                point [ x 10.0000 y 10.0000 ]
                point [ x 100.000 y 100.000 ]
                point [ x 80.0000 y 30.0000 ]
                point [ x 120.000 y 230.000 ]
                point [ x 73.0000 y 160.000 ]
              ]
            ]
          ]
        ]
    */
}
