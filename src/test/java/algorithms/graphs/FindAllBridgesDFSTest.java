package algorithms.graphs;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FindAllBridgesDFSTest extends TestCase {
    
    public FindAllBridgesDFSTest(String testName) {
        super(testName);
    }
 
    public void _test0() {
        /*
        from Figure 4 of
        lecture 4 notes of David Mount for CMSC 451 
        Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
        https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect04-edge-connectivity.pdf
        */
        
        TObjectIntMap<String> letters = new TObjectIntHashMap<String>();
        String letter;
        int i;
        int n = 10;
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[n];
        for (i = 0; i < n; ++i) {
            //A code point is a value that can be used in a coded character set.
            //A code point is a 32-bit int data type, where the lower 21 bits 
            //     represent a valid code point value and the upper 11 bits are 0.
            //A Unicode code unit is a 16-bit char value.
            
            //A is ascii 65
            //a is ascii 97
            letter = String.valueOf((char)(i + 97));
            letters.put(letter, i);
            
            g[i] = new SimpleLinkedListNode();
        }
        // add edges of vertexes
        g[letters.get("a")].insert(letters.get("e"));
        g[letters.get("a")].insert(letters.get("b"));
        g[letters.get("b")].insert(letters.get("f"));
        g[letters.get("b")].insert(letters.get("c"));
            g[letters.get("b")].insert(letters.get("a"));
        g[letters.get("c")].insert(letters.get("g"));
        g[letters.get("c")].insert(letters.get("f"));
            g[letters.get("c")].insert(letters.get("b"));
        g[letters.get("d")].insert(letters.get("h"));
        g[letters.get("d")].insert(letters.get("g"));
        g[letters.get("e")].insert(letters.get("i"));
        g[letters.get("e")].insert(letters.get("j"));
        g[letters.get("f")].insert(letters.get("b"));
        g[letters.get("f")].insert(letters.get("c"));
        g[letters.get("g")].insert(letters.get("d"));
        g[letters.get("g")].insert(letters.get("h"));
        g[letters.get("g")].insert(letters.get("c"));
        g[letters.get("h")].insert(letters.get("g"));
        g[letters.get("h")].insert(letters.get("d"));
        g[letters.get("i")].insert(letters.get("e"));
        g[letters.get("i")].insert(letters.get("j"));
        g[letters.get("j")].insert(letters.get("e"));
        g[letters.get("j")].insert(letters.get("i"));
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(0, 1));
        expected.add(new PairInt(0, 4));
        expected.add(new PairInt(2, 6));
        
        FindAllBridgesDFS dfs = new FindAllBridgesDFS(g);
        PairIntArray bridges = dfs.walk();
        assertEquals(expected.size(), bridges.getN());
        
        int u, v;
        //System.out.println("bridges:");
        for (i = 0; i < bridges.getN(); ++i) {
            u = bridges.getX(i);
            v = bridges.getY(i);
            /*System.out.printf("%d,%d, (%s,%s)\n", u,v,
               String.valueOf((char)(u + 97)),
               String.valueOf((char)(v + 97)));*/
            expected.remove(new PairInt(u, v));
        }
        assertEquals(0, expected.size());
    }

    public void test1() {
        /*
        0 -->1
        0 -->3
        1-->2
        2-->3
        2-->4

           0
           |  \
           |    1
           |     |
           |    2
           |  /  \
           3       4

        3 and 2 and 1 have only 1 incoming edge, but 2 has a predecessor that
        reaches an ancestor of 2.  1 also has a predecessor that reaches an
        ancestor of 1.

        bridges: 2:4
         */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[5];
        g[0] = new SimpleLinkedListNode();
        g[1] = new SimpleLinkedListNode();
        g[2] = new SimpleLinkedListNode();
        g[3] = new SimpleLinkedListNode();
        g[4] = new SimpleLinkedListNode();

        g[0].insert(1); g[0].insert(3);
        g[1].insert(2);
        g[2].insert(3);
        g[2].insert(4);

        // add the reverse of those edges
        g[1].insert(0);
        g[3].insert(0);
        g[2].insert(1);
        g[3].insert(2);
        g[4].insert(2);

        FindAllBridgesDFS dfs = new FindAllBridgesDFS(g);
        PairIntArray bridges = dfs.walk();
        assertTrue(bridges.getN() == 1);
        assertEquals(bridges.getX(0), 2);
        assertEquals(bridges.getY(0), 4);
    }
}
