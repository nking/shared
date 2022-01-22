package algorithms.encoding;

import algorithms.bipartite.YFastTrieWrapper;
import algorithms.heapsAndPQs.HeapNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.ArrayDeque;
import java.util.Queue;

/**
 Huffman encoding is a compression algorithm that uses entropy 
 and a prefix-free variable length code.
 <pre>
 prefix-free variable length code (a.k.a. prefix code)
    -- No codeword is a prefix of any other codeword
    -- Prefix codes are uniquely decodable
   Example: 1, 10, 110,111

   Prefix-free coding:
    -- uses a symbol tree for all symbols of message
    -- usually the most frequently used symbols are at top of tree
    -- static algorithms have 3 phases for ENCODING
        (1) read source messages and calculate symbol frequency
        (2) build symbol code tree, optimized by frequency
        (3) encode messages using symbol codes
    -- static algorithms have 3 phases for DECODING:
        (1) read the stats of the source message
            (usually the stored symbol table or other means
            to create it)
        (2) build the code tree
        (3) decode the messages
</pre>
   Huffman and Shannon-Fado attempt to assign shorter code words to more frequent
   symbols.
   
 * Huffman assigns a symbol to each
 * character based upon the frequency of that character in the text string.
 * The algorithm creates a table of frequency of characters and converts those
   to variable length binary strings.  The most frequent characters
   use the smallest binary strings.
   The compressed results are usually 20%-90% of the original size.
   
  From Cormen et al. "Introduction to algorithms":
  consider only prefix-codes, that is, codes in which no codeword is also a 
  prefix of some other codeword.
   Encoding for any binary character code is a concatenation of the codewords 
   representing each character of the file. 
   For example, with the variable-length prefix code of Figure 16.3, we code 
   the 3-character file abc as 0.101.100, where “.” denotes concatenation.
   Prefix codes are desirable because they simplify decoding. 
   Since no codeword is a prefix of any other, the codeword that begins an 
   encoded file is unambiguous. We can simply identify the initial codeword, 
   translate it back to the original character, and repeat the decoding process 
   on the remainder of the encoded file. In our example, the string 001011101 
   parses uniquely as 0 . 0 . 101 . 1101, which decodes to aabe.
 * @author nichole
 */
public class Huffman {
    
    public static class HuffmanEncoding {

        public HeapNode encodedTree;
        public TIntIntMap freqMap;
    }

    /**
     * use Huffman encoding to compress the text s. returns a tree representing
     * the compressed text and the frequency table needed to decode it. The
     * runtime complexity is O(n * log_2 log_2 (n)) where n is the length of the
     * text string c, though can be as high as O(n * log_2(n)) if Fibonacci heap
     * is used.
     *
     * @param c text to encode
     * @return returns the Huffman encoded text and the symbol code map.
     */
    public HuffmanEncoding compress(String c) {

        int n = c.length();

        // key = code point of character, value = number of times character is in c
        int sumF = 0;
        TIntIntMap f = new TIntIntHashMap();
        int cP;
        int fP;
        int i;
        for (i = 0; i < n; ++i) {
            cP = c.codePointAt(i);
            fP = 1;
            if (f.containsKey(cP)) {
                fP += f.get(cP);
            }
            f.put(cP, fP);
            sumF += fP;
        }

        // key = code point, value = binary symbol
        TIntIntMap symbolTree = buildSymbolCodeTree(f, sumF);

        // encode using symbolTree
        throw new UnsupportedOperationException("not yet complete");
    }

    protected HeapNode buildFrequencyCodeTree(TIntIntMap f, int sumF) {
        
        int n = f.size();

        // initialize heap with C set of characters which are the keys of f
        int maxQ = sumF;
        YFastTrieWrapper q = new YFastTrieWrapper(maxQ + 1);
        TIntIntIterator iter = f.iterator();
        int cP;
        int fP;
        // or PriorityQueue<Integer> q = new PriorityQueue<>(n);
        HeapNode node;
        while (iter.hasNext()) {
            iter.advance();;
            cP = iter.key();
            fP = iter.value();
            node = new HeapNode(fP);
            node.setData(cP);
            q.insert(node);
        }

        //build a tree from the bottom up
        HeapNode nodeZ;
        int i;
        for (i = 1; i < n; ++i) {
            nodeZ = new HeapNode();
            nodeZ.setLeft(q.extractMin());
            nodeZ.setRight(q.extractMin());
            nodeZ.setKey(nodeZ.getLeft().getKey() + nodeZ.getRight().getKey());
            q.insert(nodeZ);
        }

        assert (q.getNumberOfNodes() == 1);

        HeapNode t = q.extractMin();

        return t;
    }

    protected TIntIntMap buildSymbolCodeTree(TIntIntMap f, int sumF) {

        HeapNode t = buildFrequencyCodeTree(f, sumF);

        TIntIntMap symbolCodeMap = new TIntIntHashMap();

        // convert the tree into variable-length binary prefix-code.
        // traversal is level-order.
        HeapNode node = t;
        // the node key here is changed from frequency to binary code
        t.setKey(0);// key=freq->code, value=cP
        Queue<HeapNode> queue = new ArrayDeque<>();
        long key;
        while (node != null) {
            if (node.getLeft() == null && node.getRight() == null) {
                // this is a leaf, so is one of the characters
                symbolCodeMap.put((Integer) node.getData(), (int) node.getKey());
            } else {
                if (node.getLeft() != null) {
                    //change left key: append 0 to node key
                    key = node.getKey() << 1;
                    node.getLeft().setKey(key);
                    queue.add(node.getLeft());
                }
                if (node.getRight() != null) {
                    //change right key: append 1 to node key
                    key = node.getKey() << 1;
                    key |= 1;
                    node.getRight().setKey(key);
                    queue.add(node.getRight());
                }
            }
            node = queue.poll(); // returns null if empty
        }

        return symbolCodeMap;
    }

    /**
     * decode the Huffman encoded tree. runtime complexity is O(N) where N is
     * the number of characters in the original text.
     *
     * @param h the Huffman encoded text and the symbol code map
     * @return returns the decoded text.
     */
    public String decompress(HuffmanEncoding h) {
        throw new UnsupportedOperationException("not yet complete");
    }
}
