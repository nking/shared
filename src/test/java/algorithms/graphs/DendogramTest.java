package algorithms.graphs;

import algorithms.misc.MinMaxPeakFinder;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
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
        TFloatList qs = new TFloatArrayList();
        
        Modularity m = new Modularity();
        
        double[][] e;
        int i, j, k;
        double q;
        for (Dendogram.DendogramLayer layer : layers) {
            
            q = m.girvanNewman2002(layer, adjList, nEdges);
            k = layer.nComponents;
            
            ks.add(k);
            qs.add((float)q);
        }
        
        for (i = 0; i < ks.size(); ++i) {
            System.out.printf("k=%d q=%.3f%n", ks.get(i), qs.get(i));
        }
        
        MinMaxPeakFinder mpf = new MinMaxPeakFinder();
        int[] maxIdxs = mpf.findPeaks(qs.toArray(), 0.f, 2.5f);
        for (i = 0; i < maxIdxs.length; ++i) {
            int idx = maxIdxs[i];
            System.out.printf("found peak for q=%.3f k=%d%n", qs.get(idx), ks.get(idx));
        }
        
    }

}
