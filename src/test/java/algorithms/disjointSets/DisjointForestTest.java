package algorithms.disjointSets;

import algorithms.disjointSets.DisjointForest.RootedTreeDisjointSet;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.set.TIntSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DisjointForestTest extends TestCase {
    
    public DisjointForestTest(String testName) {
        super(testName);
    }

    public void test0() {
        // from Fig 21.4 of Corment et al. Intro to Alg.
        
        // b, c, e, h  is one tree
        DisjointSet2Node<String> b = new DisjointSet2Node<String>("b");
        DisjointSet2Node<String> c = new DisjointSet2Node<String>("c");
        DisjointSet2Node<String> e = new DisjointSet2Node<String>("e");
        DisjointSet2Node<String> h = new DisjointSet2Node<String>("h");
        
        
        DisjointForest<String> forest = new DisjointForest<>();
        
        DisjointForest.RootedTreeDisjointSet<String> treeB = forest.makeSet(b);
        DisjointForest.RootedTreeDisjointSet<String> treeC = forest.makeSet(c);
        DisjointForest.RootedTreeDisjointSet<String> treeE = forest.makeSet(e);
        DisjointForest.RootedTreeDisjointSet<String> treeH = forest.makeSet(h);
        assertEquals(4, forest.getTrees().size());
        
        System.out.println("\nforest(b,c,e,h) before any unions:\n" + forest.toString());
        
        DisjointSet2Node<String> bh = forest.union(b, h);
        System.out.println("\nforest(b,c,e,h) after any union b,h:\n" + forest.toString());
        assertEquals(3, forest.getTrees().size());
        
        DisjointSet2Node<String> ec = forest.union(e, c);
        System.out.println("\nforest(b,c,e,h) after any union e,c:\n" + forest.toString());
        assertEquals(2, forest.getTrees().size());
        
        DisjointSet2Node<String> bhec = forest.union(b, e);
        assertEquals(1, forest.getTrees().size());
        System.out.println("\nforest(b,c,e,h) after any union all:\n" + forest.toString());
        
        // findSet compresses c's path to point to c's earliest ancestor parent
        DisjointSet2Node<String> parentPOfC = forest.findSet(c);
        System.out.println("\nforest(b,c,e,h) after findSet(c):\n" + forest.toString());
        assertEquals(1, forest.getTrees().size());
        
        Iterator<Entry<DisjointSet2Node<String>, DisjointForest.RootedTreeDisjointSet<String>>> iter 
            = forest.getTrees().entrySet().iterator();
        while (iter.hasNext()) {
            Entry<DisjointSet2Node<String>, DisjointForest.RootedTreeDisjointSet<String>>
                entry = iter.next();
            assertEquals(parentPOfC, entry.getKey());
            assertEquals(parentPOfC, entry.getValue().parent);
        }
    }
    
    public void test1() {
        /*
        public static Map<DisjointSet2Node<Integer>, RootedTreeDisjointSet<Integer>> 
        connectedComponents(SimpleLinkedListNode[] adjList) {
        */
        int a = 0; int b = 1; int c = 2; int d = 3;int e = 4; 
        int f = 5; int g = 6; int h = 7;
        SimpleLinkedListNode[] adjList = new SimpleLinkedListNode[8];
        for (int i = 0; i < 8; ++i) {
            adjList[i] = new SimpleLinkedListNode();
        }
        int nUnassigned = 1;
        
        adjList[c].insert(h); adjList[c].insert(e);
        adjList[h].insert(b);
        
        adjList[f].insert(d);
        adjList[d].insert(g);
        
        List<TIntSet> comp = DisjointForest.connectedComponents(adjList);
        
        assertEquals(2 + nUnassigned, comp.size());
    }
}
