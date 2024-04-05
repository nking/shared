package algorithms.encoding;

import algorithms.TreeTraversal;
import algorithms.VeryLongBitString;
import algorithms.bipartite.YFastTrieWrapper;
import algorithms.heapsAndPQs.HeapNode;
import algorithms.misc.MiscMath0;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.math.BigInteger;
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
   
  From Cormen, Leiserson, Rivest, and Stein "Introduction to algorithms":
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

    /**
     *
     */
    public static class HuffmanEncoding {
        
        /**
         * The symbol code tree:
         * the values in the leaf nodes as datum are the code point for decoded letters.
         * To use the tree with the encoded bit-string: start from the root
         * node, and read the first bit in the bit-string.  A '0' in the bit-string
         * means follow the left node while a '1' in the bit-string means
         * follow the right node.  When one reaches a leaf, that is the
         * letter that has just been decoded, so write it.  Continue reading the bit-string
         * but begin the traversal from the root node again.
         */
        public HeapNode symbolTree = null;
        
        /**
         * a map with key=code point of a letter in the alphabet to encode,
         * and the value for the key is the number of times the letter
         * occurs in the text to be encoded.
         * This can be used to build the symbol code tree,
         */
        public TIntIntMap freqMap = null;
        
        /**
         * the encoded text.  when reading, start from bit 0.
         */
        public VeryLongBitString encoded = null;
        
        /**
         *
         */
        public void print() {
            System.out.println("freqMap:");
            if (freqMap != null) {
                printMap(freqMap);
            }

            System.out.println("symbolTree:");
            if (symbolTree != null) {
                TreeTraversal.printLevelOrder(symbolTree);
            }
            
            System.out.println("encoded bitstring:");
            if (encoded != null) {
                System.out.println(encoded.toString());
            }
        }
    }
    
    /**
     *
     */
    protected static class EncodingSymbols {
        
        /**
         * The symbol code tree:
         * the values in the leaf nodes as datum are the code point for original uncoded letters
         * while the keys are the bit-strings for the encodings.
         * Note that the nodes which are not leaves (that is not encoded alphabet)
         * will have null for their data property.
         */
        HeapNode symbolTree = null;
        
        /**
         * keys are the code-points for the original un-coded text.
         * values for the keys are the bit-strings for the encodings.
         */
        TIntIntMap codeSymbolMap = null;
        
        /**
         *
         */
        public void print() {
            System.out.println("codeSymbolMap:");
            if (codeSymbolMap != null) {
                printMap(codeSymbolMap);
            }

            System.out.println("symbolTree:");
            if (symbolTree != null) {
                TreeTraversal.printLevelOrder(symbolTree);
            }
        }
    }

    /**
     * use Huffman encoding to compress the text s. returns a tree representing
     * the compressed text and the frequency table needed to decode it. The
     * runtime complexity is O(n * log_2 log_2 (n)) where n is the length of the
     * text string c, though can be as high as O(n * log_2(n)) if Fibonacci heap
     * is used.
     *
     @param c text to encode
     @return returns the Huffman encoded text and the symbol code map.
     */
    public HuffmanEncoding compress(String c) {

        int n = c.length();

        // key = code point of character, value = number of times character is in c
        TIntIntMap f = new TIntIntHashMap();
        int sumF = 0;
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

        EncodingSymbols symbols = buildSymbolCodeTree(f, sumF);
        
        // encode c         
        VeryLongBitString encoded = encode(c, symbols.codeSymbolMap);
        
        HuffmanEncoding h = new HuffmanEncoding();
        h.symbolTree = symbols.symbolTree;
        h.freqMap = f;
        h.encoded = encoded;
        
        return h;
    }

    /**
     @param f map with key=code-point, value=frequency of code-point in 
     * original text.
     @param sumF the sum of all frequencies in f.  it's used to assign
     * a maximum bit-length needed by an internal YFastTrie data structure.
     @return 
     */
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
            iter.advance();
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
        
        t = repairToUniqueIfNeeded(t);
        
        return t;
    }
    
    /**
     * 
     @param f map with key=code-point, value=frequency of code-point in 
     * original text.
     @param sumF the sum of all frequencies in f.  it's used to assign
     * a maximum bit-length needed by an internal YFastTrie data structure.
     @return 
     */
    protected EncodingSymbols buildSymbolCodeTree(TIntIntMap f, int sumF) {
        
        //System.out.println("freqMap:");
        //printMap(f);   
        
        // key = combined frequency, data=code point.
        HeapNode t = buildFrequencyCodeTree(f, sumF);
        
        
        //System.out.println("freq code tree:");
        //TreeTraversal.printLevelOrder(t);
        

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
        
        EncodingSymbols eS = new EncodingSymbols();
        eS.symbolTree = t;
        eS.codeSymbolMap = symbolCodeMap;
        
        //eS.print();

        return eS;
    }
    
    private long count(String c, TIntIntMap codeSymbolMap) {
        long nBits = 0;
        int cP;
        int s;
        for (int i = 0; i < c.length(); ++i) {
            cP = c.codePointAt(i);
            s = codeSymbolMap.get(cP);
            nBits += MiscMath0.numberOfBits(s);
        }
        return nBits;
    }

    /**
     * encode the text c, given the code symbol map.
     @param c uncoded text to be encoded.
     @param codeSymbolMap map with key = alphabet code point, value = 
     * bit-string symbol.
     @return returns the encoding of text 'c', stored as a VeryLongBitString.
     */
    VeryLongBitString encode(String c, TIntIntMap codeSymbolMap) {
        
        // count the number of bits the encoded string will be.
        // this is needed to construct the VeryLongBitString.
        long nE = count(c, codeSymbolMap);
        
        //NOTE: could use java's BigInteger, but that creates a new BigInteger 
        //      for every modifying operation such as setBit().
        VeryLongBitString encoded = new VeryLongBitString(nE);

        int cP;
        int s;
        long b = 0;
        long b2;
        int nB;
        for (int i = 0; i < c.length(); ++i) {
            cP = c.codePointAt(i);
            s = codeSymbolMap.get(cP);
            
            nB = MiscMath0.numberOfBits(s);

            //TODO: revisit one day
            
            // set bits, reversing them from b to b+nB for the bit-string tree path
            // _ _ _ _ _ _  b=0, nB=3   7 6 5 4 3 2 1 0
            //              s=100: write s[0] to  c[(b+nB-1)]
            b2 = b + nB - 1;
            while (nB > 0) {
                if ((s & 1) == 1) {
                    encoded.setBit(b2);
                }
                s >>= 1;
                ++b;
                --nB;
                --b2;
            }
        }
        
        return encoded;
    }
    
    /**
     * decode the Huffman encoded tree. runtime complexity is O(N) where N is
     * the number of characters in the original text.
     *
     @param symbolTree the symbol code map
     @param encoded text
     @return returns the decoded text.
     */
    public String decompress(HeapNode symbolTree, VeryLongBitString encoded) {
        
        StringBuilder sb = new StringBuilder();
        
        HeapNode node = symbolTree;
             
        if (node == null) {
            return "";
        }
        
        long n = encoded.getInstantiatedBitSize();
        long i;
        int cP;
        for (i = 0; i < n; ++i) {
            if (node.getLeft() == null && node.getRight() == null) {
                
                // this is a leaf node, so we have a letter in the decoded alphabet now
                cP = (Integer)node.getData();
                
                sb.append(Character.toChars(cP));
                node = symbolTree;
            }
            
            if (encoded.isSet(i)) {
                node = node.getRight();
            } else {
                node = node.getLeft();
            }
        }
        
        if (node != null && node.getData() != null) {
            cP = (Integer)node.getData();
            sb.append(Character.toChars(cP));
        }
        
        return sb.toString();
    }
    
    private HeapNode repairToUniqueIfNeeded(HeapNode t) {
        
        //System.out.println("\nrepairToUniqueIfNeeded:");
        //    TreeTraversal.printLevelOrder(t);
            
        // there is a possible violation of unique prefix code for
        // a path from root to leaf that is all left nodes and with a path length longer than 1.
        // can remove the "all-left" path leaf
        // add a new root node to top of the tree t and assign the existing tree to its right
        // assign the leaf to the new root's left 
        
        if (!(t.getLeft() != null && t.getLeft().getLeft() != null)) {
            return t;
        }
       
        HeapNode p = t;
        HeapNode leaf = p.getLeft();
                
        while (leaf.getLeft() != null) {
            p = leaf;
            leaf = leaf.getLeft();
        }
        
        //System.out.printf("remove leaf=(%d,%d), p=(%d,%d)\n", 
        //    leaf.getKey(), (Integer)leaf.getData(), p.getKey(), (Integer)p.getData());
        p.setLeft(null);
        
        HeapNode r = new HeapNode(t.getKey() - 1);
        r.setRight(t);
        r.setLeft(leaf);
        
        //System.out.println("repaired tree:");
        //TreeTraversal.printLevelOrder(r);
        
        return r;
    }
    
    private static void printMap(TIntIntMap f) {
        TIntIntIterator iter = f.iterator();
        while (iter.hasNext()) {
            iter.advance();
            System.out.printf("key=%d, v=%d\n", iter.key(), iter.value());
        }
    }

}
