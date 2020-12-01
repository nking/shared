package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
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
        
        double nEdges = g.edgeWeightMap.size();
        
        // store q and k to look for peak in q
        TIntList ks = new TIntArrayList();
        TDoubleList qs = new TDoubleArrayList();
        
        double[][] e;
        int i, j, k, c1, c2;
        double trE, q;
        for (Dendogram.DendogramLayer layer : layers) {
            
            k = layer.nComponents;
            
            //create matrix e where e[i][j] is fraction of edges from i to jover all edges
            e = new double[k][k];
            for (i = 0; i < k; ++i) {
                e[i] = new double[k];
            }
            
            // count number of edges between communities using direction of edges
            for (i = 0; i < layer.vertexComponents.length; ++i) {
                c1 = layer.vertexComponents[i];
                if (c1 == -1) {
                    continue;
                }
                
                SimpleLinkedListNode jNode = adjList[i];
                while (jNode != null && jNode.getKey() != -1) {
                    j = jNode.getKey();
                    c2 = layer.vertexComponents[j];
                    if (c2 == -1) {
                        continue;
                    }
                    e[c1][c2]++;
                    jNode = jNode.getNext();
                }
            }
            
            // divide matrix e by number of edges
            for (i = 0; i < k; ++i) {
                for (j = 0; j < k; ++j) {
                    e[i][j] /= nEdges;
                }
            }
            
            //Trace of e is the fraction of edges in the network that connect 
            //   vertices in the same community.  this is high for good component divisions.
            trE = 0;
            for (i = 0; i < k; ++i) {
                trE += e[i][i];
            }
            
            // Q = trE - || e^2 ||
            q = 0;
            for (i = 0; i < k; ++i) {
                for (j = 0; j < k; ++j) {
                    q += (e[i][j] * e[i][j]);
                }
            }
            q = trE - q;
            
            ks.add(k);
            qs.add(q);
        }
        
        for (i = 0; i < ks.size(); ++i) {
            System.out.printf("k=%d q=%.3f\n", ks.get(i), qs.get(i));
        }
    }

}
