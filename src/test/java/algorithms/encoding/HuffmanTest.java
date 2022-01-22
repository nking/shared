package algorithms.encoding;

import algorithms.heapsAndPQs.HeapNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HuffmanTest extends TestCase {
    
    public HuffmanTest(String testName) {
        super(testName);
    }

    public void testCompress() {
        
    }

    /**
     * Test of buildFrequencyCodeTree method, of class Huffman.
     */
    public void testBuildFrequencyCodeTree() {
        
        // Fig 16.4 from Cormen et al. "Introduction to Algorithms"
        
        System.out.println("buildFrequencyCodeTree");
        
        TIntIntMap f = new TIntIntHashMap();
        f.put("a".codePointAt(0), 45);
        f.put("b".codePointAt(0), 13);
        f.put("c".codePointAt(0), 12);
        f.put("d".codePointAt(0), 16);
        f.put("e".codePointAt(0), 9);
        f.put("f".codePointAt(0), 5);
        
        int sumF = 45 + 13 + 12 + 16 + 9 + 5;
        
        Huffman h = new Huffman();
        HeapNode t = h.buildFrequencyCodeTree(f, sumF);
        
        assertEquals(100, t.getKey());
        assertEquals(45, t.getLeft().getKey());
        assertEquals(55, t.getRight().getKey());
        assertEquals(25, t.getRight().getLeft().getKey());
        assertEquals(12, t.getRight().getLeft().getLeft().getKey());
        
        assertEquals(30, t.getRight().getRight().getKey());
        assertEquals(16, t.getRight().getRight().getRight().getKey());
    
        TIntIntMap s = h.buildSymbolCodeTree(f, sumF);
        
        assertEquals(6, s.size());
        assertEquals(0, s.get("a".codePointAt(0)));
        assertEquals(5, s.get("b".codePointAt(0)));
        assertEquals(4, s.get("c".codePointAt(0)));
        assertEquals(7, s.get("d".codePointAt(0)));
        assertEquals(13, s.get("e".codePointAt(0)));
        assertEquals(12, s.get("f".codePointAt(0)));
    }

    /**
     * Test of decompress method, of class Huffman.
     */
    public void testDecompress() {
        System.out.println("decompress");
    }
    
}
