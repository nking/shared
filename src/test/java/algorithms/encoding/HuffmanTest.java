package algorithms.encoding;

import algorithms.VeryLongBitString;
import algorithms.encoding.Huffman.HuffmanEncoding;
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

    /**
     * Test of buildFrequencyCodeTree method, of class Huffman.
     */
    public void testBuildFrequencyCodeTree() {
        
        // Fig 16.4 from Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"
        
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
    
        Huffman.EncodingSymbols symbols = h.buildSymbolCodeTree(f, sumF);
        TIntIntMap s = symbols.codeSymbolMap;
        
        assertEquals(6, s.size());
        assertEquals(0, s.get("a".codePointAt(0))); //0
        assertEquals(5, s.get("b".codePointAt(0))); //101
        assertEquals(4, s.get("c".codePointAt(0))); //100
        assertEquals(7, s.get("d".codePointAt(0))); //111
        assertEquals(13, s.get("e".codePointAt(0)));//1101
        assertEquals(12, s.get("f".codePointAt(0)));//1100
        
        assertEquals("a".codePointAt(0), ((Integer)symbols.symbolTree.getLeft().getData()).intValue());
        
        assertEquals("b".codePointAt(0), 
            ((Integer)symbols.symbolTree.getRight().getLeft().getRight().getData()).intValue());
        
        assertEquals("c".codePointAt(0), 
            ((Integer)symbols.symbolTree.getRight().getLeft().getLeft().getData()).intValue());
        
        assertEquals("d".codePointAt(0), 
            ((Integer)symbols.symbolTree.getRight().getRight().getRight().getData()).intValue());
        
        assertEquals("e".codePointAt(0), 
            ((Integer)symbols.symbolTree.getRight().getRight().getLeft().getRight().getData()).intValue());
        
        assertEquals("f".codePointAt(0), 
            ((Integer)symbols.symbolTree.getRight().getRight().getLeft().getLeft().getData()).intValue());
        
        //"a"//0
        //"b"//101
        //"c"//100
        //"d"//111
        //"e"//1101
        //"f"//1100
        
        String uncoded = "cbbfaed";
        
        String e = "111101100011101101001";
        
        VeryLongBitString encoded = h.encode(uncoded, s);        
        System.out.printf("encoded= %s\n", encoded);
        System.out.printf("expected=%s\n", e);
        
        for (int i = 0; i < e.length(); ++i) {
            if (e.charAt(e.length() - 1 - i) == '1') {
                assertTrue(encoded.isSet(i));
            } else {
                assertTrue(encoded.isNotSet(i));
            }
        }
        System.out.println();
        
        String decoded = h.decompress(symbols.symbolTree, encoded);
        
        System.out.println(decoded);
        assertTrue(uncoded.equals(decoded));
        
        HuffmanEncoding he = h.compress(uncoded);
        System.out.printf("encoded=%s;  path=%s\n", he.encoded.toString(),
            new StringBuilder(he.encoded.toString()).reverse().toString());
        decoded = h.decompress(he.symbolTree, he.encoded);
        System.out.println(" compress, decompress =" + decoded);
        assertTrue(uncoded.equals(decoded));
    }
    
    public void testCompressDecompress() {
        
        // from wikipedia https://en.wikipedia.org/wiki/Photophone
        String uncoded = "The photophone is a telecommunications device that "
            + "allows transmission of speech on a beam of light. It was "
            + "invented jointly by Alexander Graham Bell and his assistant "
            + "Charles Sumner Tainter on February 19, 1880, at Bell's "
            + "laboratory at 1325 L Street in Washington, D.C.";
        
        Huffman h = new Huffman();
        
        HuffmanEncoding he = h.compress(uncoded);
        
        String decoded = h.decompress(he.symbolTree, he.encoded);
        
        System.out.println("decoded=" + decoded);
        
        assertTrue(uncoded.equals(decoded));
        
    }

}
