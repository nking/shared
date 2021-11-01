package algorithms.tsp;

import algorithms.SubsetChooser;
import algorithms.misc.MiscMath0;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Arrays;

/**
 * * <pre>
 * 
 * Keeping notes here on looking at bitstring patterns to find ways to 
 * condense the sums of 3-paths used in the subproblems
 * for a dynamic solution to TSP as an exercise in recursion and
 * dynamic programming (not a preferred impl due to the large runtime complexity).
 * 

    let start node = 0.  n=5
    The total number of permutations of a path with fixed start node = 0 and n=5
    is n!/n = 24.

    for n=127, nPerm = 2.4e+211, so if wanted to store each permutation, indexed,
    would need to use java's BigInteger as a bitstring or this project's VeryLongBitString.java

    to re-use sub-problems, can determine how many sets of paths of 3 nodes are in
    the set of n-1 which is (n-1)/3.
    (sets of size 2 are given by the distance matrix.)

    Store all 3 node combinations, excluding the number 0, for reuse.
    The number of k-permutations without excluding a node is: C(n,k) = n!/(k!*(n-k)!).
    For example, 4 nodes, k-permutations size=3:  C(n,k) = n!/(k!*(n-k)!) = 4
      7 (    111)
     11 (   1011)
     13 (   1101)
     14 (   1110)
    Excluding the permutations with bit=0 set: use n2=n-1, C(n-1,k) = 1
     14 (   1110)

    So, one can generate the k=3 permutations of subset of size n, and exclude the 1st bit
    by generating for n2=n-1, and the left shift by 1 the result.
    (can use SubsetChooser.java for the permutations)

    To further look at bit patterns to see if memoization can be used, but with
    condensed intermediate steps to use less storage, one more example:

    For n=10, nPerm=84:
    as the 3 sequence bitstrings are generated, the inverse bitstring w/ 0 cleared
    has the set of nodes this bitstring can be added to
         bitstring    inverse, excl 0
    14 (0000001110)  (1111110000)
    22 (0000010110)  (1111101000) perm for n3=125 can be generated, exclude bit 1, then left shift by 3 to get paths to add this to.
    26 (0000011010)  (1111100100) perm for n3=249 can be generated, exclude bits 1,2, then left shift by 2 to get paths to add this to.
    28 (0000011100)  (1111100000)
    38 (0000100110)  (1111011000)
    42 (0000101010)  (1111010100)
    44 (0000101100)  (1111010000)
    50 (0000110010)  (1111001100)
    52 (0000110100)  (1111001000)
    56 (0000111000)  (1111000100)
    70 (0001000110)  (1110111000)
    74 (0001001010)  (1110110100)
    ...
    n=10, k=3 => count = 84

    Once all of the 3-node path subsets are generated and stored,
    one can consider again that each of the complete permutations of paths
    is composed of disjoint combinations of (n-1)/3 of the 3-node path subsets,
    (plus up to 2 nodes.  e.g. if n = 12, we have (n-1)/3 3-node path subsets
     + 1 fixed node + 2 free nodes that belong in the permuation).

    For now, will consider only the cases for n = 1 + a multiple of 3 to look at
    ways to condense the problem.

    for the n=10 exaple, we have 84 k=3 sequences and have the unset bits, excluding bit 0

    note: as each complete path of the (n-1)/3 sections is totaled, compare it to minTotal.

    at this point, we have generated and stored all 3-node path sums, excluding the start node.

    we could proceeed in many ways:
      (1) one way would be to use recursion all the way to the end of a path, and add the
          3-node path sums to each stage, but not save the intermediate stages.
          the memory requirements would be more manageable.
          The sub-problem re-use would be only of the 3-node paths though, and could not re-use
          already summed partial paths longer than the 3-nodes.
          Note that this pattern is just one level from a brute force approach.
          (and in that case, one might want to consider designing an algorithm composed of
          as much memoization as possible for the system the code is running on while completing
          the rest of the calculations with recursive brute force approach.  e.g.
          enough memory to store the 3-node permutations?  6-node permutations?  9-node permutations?
          then continue each path to its end with recursion and compare its total 
          with mintotal and discard it if not mincost).
      (2) another way to proceed would be to use recursion all the way to the end of a path,
          and to save all intermediate sums for re-use.
          The sub-problem re-use is then truly the dynamic approach, but requires more memory
          for storage.

          The dynamic approach needs to store the partial sums as the path lengthens
             (without the start node).
          datastructure requirements:
             - need set(key, value), getValue(key), and contains(key), and prefer O(1) or O(small const) for each.
             - a key for one partial sum would have to hold information on the nodes in
               the partial path and their order.
               could consider many ways to do that:
                 The number of bits needed to represent the number of a node would be math.log(n-1)/math.log(2).
                 Let w = math.log(n-1)/math.log(2).  then a key as a bitstring representing node and order
                 would be a variable bit length from 3 nodes = 3*w bits in length up to (n-1)*w bits in length.
                 Ror example, a 3-node path sum from nodes 2, 3, 1 would be
                      2       3       1
                   ------  ------  ------
                   w bits  w bits  w bits

                 Compression options for the key could be considereed too.

               Note that the start node isn't part of those keys as it will be added into the path at
               evaulation time.
             - A worse case number of keys to store is the permutation of all but one node: n!/n.
             - A long bitstring is limited to 63 bits.
             - A java array length is limited to signed integer length, 1 &lt&lt 31 -1.
             - so to fit more than 1 &lt&lt 31-1 permutations would need to use the bitstring key
               as a concatenation of more than one path bitstring.
               the java BigInteger or the VeryLongBitString could do so.
               - also, the values for those keys would need to be referenced in an O(1) manner also if possible.
                 one could use a 2 dimensional array where the 2nd dimension
                 holds the individual values for the keys in the concatenated key.
             - would need a wrapper class to manage the integrity of the array of BigInteger or
               VeryLongBitString and the 2 dimensional array of path sums.
             - what amount of concatenation of keys would improve the number of storable
               path sums significantly?

       Will continue from here having chosen (2) to complete the dynamic approach.
       Will also assume that the datastructure for storage exists and will call it memo for now.
       
       considering recursion patterns:
       
           calculateAndStore3NodePaths();
           int nNodesRemaining = (n-1) - 3;
           r3(bitstring, sum, nNodesSet);
           return min path(s) and the min cost

           private void r3(bitstring, sum, int nNodesRemaining) {
               //if possible, use tail recursion in design...
               //    (best for C++, java doesn't use tail recursion for method frames)
               if (nNodesRemaining == 0) {
                   add start and end node costs
                   compare to min and if smaller or same, store min cost and path
                   (possibly would like to store all min paths).
                   note that there could be user option to set tolerance in comparison of cost being the same.
                   return;
               }
               if (nNodesRemaining <= k) {
                   // assert k = 3.  if that changes from tinkering, the assert will alarm that this section will fail.
                   calculate the 1 or 2 permutations of those remaining 1 or 2 nodes
                   calc the concatenated bit strings and the the path sums
                   store those in memo
                   invoke r3 for the concatenated bitstrings each
                   return;
               }
               int nBitsSet = (n-1) - nNodesRemaining;
               subsetchooser = new...(nNodesRemaining, k)
               while (true) {
                   s = subsetchooser.next();
                   if (s == -1){break;}
                   int[] si = tranform s to bitstring unset indexes.  should be 3 bits set
                   bitstring2 = concatenate(bitstring, nBitsSet, si);
                   sum2 = sum + memo.get(si); <== edit here.  need a method to convert si to a memo bitstring
                   memo.set(bitstring2, sum2);
                   r3(bitstring2, sum2, nNodesRemaining - 3);
               }
           }
           
rewrite in iterative form:

           calculateAndStore3NodePaths();
           int nNodesRemaining = (n-1) - 3;
           min = Long.POSITIVE_INFINITY;
           minPaths = null;
           for each bitstring in the 3 node paths just calculated {
           
               stack = new stack;// specialized stack designed for pairs of bitstring keys and sum and nNodesRemaining,
                                 // that uses same pattern as memo to compress keys to allow more items to be stored
                                 
               stack.add(bitstring, sum, nNodesRemaining);

               while (!stack.isEmpty()) {

                   bitstring2, sum2, nNodesRemaining2 = stack.pop();

                   if (nNodesRemaining2 == 0) { 
                       //add start and end node costs
                       //compare to min and if smaller or same, store min cost and path
                       compareToMin(bitstring3, sum3);
                       continue;
                   }
                   if (nNodesRemaining2 <= k) {
                       // assert k = 3.  if that changes from tinkering, the assert will alarm that this section will fail.
                       calculate the 1 or 2 permutations of those remaining 1 or 2 nodes
                       calc the concatenated bit strings and the the path sums
                       store those in memo
                       compare each to min
                       continue;
                   }

                   int nBitsSet = (n-1) - nNodesRemaining2;
                   subsetchooser = new...(nNodesRemaining2, k)

                   while (true) {
                       s = subsetchooser.next();
                       if (s == -1){
                           break;
                       }
                       int[] si = tranform s to bitstring2 unset indexes.  should be 3 bits set
                       bitstring3 = concatenate(bitstring2, nBitsSet, si);
                       sum3 = sum2 + memo.get(si); <== edit here.  need a method to convert si to a memo bitstring
                       memo.set(bitstring3, sum3);

                       stack.add(bitstring3, sum3, nNorderRemaining - 3);
                   }
               }
           }
    
 * </pre>
 * 
 * @author nichole
 */
public class TSPDynamic {
   
    private long minCost = Long.MAX_VALUE;
    private TLongList minPath = new TLongArrayList();
    private final int startNode = 0;
    private final int[][] dist;
    private final long[][] memo;
    private final long sentinel = Long.MAX_VALUE;
    private final long totalNPerm;
    private final long totalNSubSeq;
    private final int w; // number of bits a city takes in a path where path is a bitstring of type long
    
    public TSPDynamic(int[][] dist) {
        
        this.dist = dist;
        int n = dist.length;
        long nPerm = MiscMath0.factorial(n); // max for n=13 for limit of array length
        totalNPerm = nPerm/n;
        totalNSubSeq = countTotalNumSubSeqInvocations(n); // max for n=507 for limit of array length
        
        if (totalNSubSeq > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("this class can solve for 13 cities at most."
                + " one could design a datastructure to hold more than (1<<31)-1 items, "
                + " but the algorithms in this class have exponential runtime complexity,"
                + " so it would be better to choose a good approximate TSP.");
        }
        n = (int)totalNPerm;
                
        w = (int)(Math.ceil(Math.log(dist.length-1)/Math.log(2)));

        memo = new long[n][];
        for (int i = 0; i < n; ++i) {
            memo[i] = new long[n];
            Arrays.fill(memo[i], sentinel);
        }
        
        /*
        //find max n cities for limit of array length
        int nc = 100;
        long wc, nb;
        while (true) {
            long c = countTotalNumSubSeq(nc);
            if (c > Integer.MAX_VALUE) {
                break;
            }
            wc = (long)Math.ceil(Math.log(nc-1)/Math.log(2)); // 1 city encoding in bits
            nb = (Long.MAX_VALUE/wc); // the number of cities one could fit into a long if using wc
            System.out.printf("nc=%d, ns=%d wc=%d nb=%d ns*1e-9=%e\n", nc, c, wc, nb, c*1.e-9);
            nc = (int)Math.round(nc*1.5);
        }
        */
    }
    
    public void test0() {
        
        //nP = C(n-1,k) for 1 fixed node
        int n = 10;
        int k = 3;
        
        long nPerm = MiscMath0.factorial(n);
        nPerm = nPerm/n;
        long nTotalPerm = countTotalNumSubSeqInvocations(n);
        System.out.printf("nPerm=%d, nTotalPerm=%d\n", nPerm, nTotalPerm);
        
        
        SubsetChooser chooser = new SubsetChooser(n-1, k);
        int c = 0;
        long s, sInv;
        StringBuilder sb;
        //https://www.geeksforgeeks.org/program-to-invert-bits-of-a-number-efficiently/
        long all1s = 1 << (n-1);
        all1s = all1s|all1s-1;
        int nSetBits;
        long tZeros;
        while (true) {
            s = chooser.getNextSubset64Bitstring();
            if (s == -1) {
                break;
            }
            
            s <<= 1;
            sb = new StringBuilder(Long.toBinaryString(s));
            
            while (sb.length() < n) {
                sb = sb.insert(0, "0");
            }
            
            sInv = s^all1s;
            sInv &= ~1; // clear bit 0
            
            nSetBits = Long.bitCount(sInv);
            tZeros = Long.numberOfTrailingZeros(sInv);
            
            System.out.printf("%d (%7s)  (%7s) nSetBits=%d tz=%d\n", 
                s, sb, Long.toBinaryString(sInv), nSetBits, tZeros);
            
            // this permutation of nSetBits can be used for all 84 of the first
            //   3-node paths.
            //   each use of the permutation can be left shifted by tZeros
            //   then assigned positions relative to the set bits in sInv
            //   to get the remaining permutations for this specific si
            SubsetChooser chooseri = new SubsetChooser(nSetBits,  k);
            long si;
            while (true) {
                si = chooseri.getNextSubset64Bitstring();
                if (si == -1) {
                    break;
                }
                long s2 = si << tZeros;
                
                // the nSetBits in si are to be shifted by tZeros
                // then converted to the si set bit positions
                // 
                
            }
            
            c++;
        }
        System.out.printf("n=%d count=%d\n", n, c);
    }
    
    // total number of subchooser invocations.  max n = 507 for count < Integer.MAX_VALUE
    protected long countTotalNumSubSeqInvocations(int n) {
        int k=3, n2=n;
        long c = 0;
        while (n2 > k) {
            c += MiscMath0.computeNDivKTimesNMinusK(n2, k);
            n2 -= k;
        }
        return c;
    }
    
    /**
     * given a pathNodeNumber such as 0 being the first node in the path bit-string,
     * uses path node bit-string length w to read the set bits and 
     * return the base 10 number (which can be used with the distance matrix 
     * or in writing the final solution of ordered nodes).
     * @param path
     * @param pathNodeNumber number of the node within the path.  NOTE: this excludes
     * start node 0 and end node 0 (pathNodeNumber=0 corresponds to the
     * 2nd node in the final solution for the completed path for the algorithm instance).
     * @return 
     */
    protected int getBase10NodeIndex(final long pathNodeNumber, final long path) {
        // read bits pathNodeNumber*w to pathNodeNumber*w + w
        long sum = 0;
        long b = pathNodeNumber*w;
        long end = b + w;
        int s = 0;
        for (b = pathNodeNumber*w, s=0; b < end; ++b, s++) {
            if ((path & (1L << b)) != 0) { // test bit b in path
                sum += (1L << s);
            }
        }
        return (int)sum;
    }
    
    /**
     * write the base 10 indexes s into a bit-string in the encoding used by the
     * memo.
     * @param s base 10 node indexes in the order to set into the bit-string path.
     * NOTE that s should not contain the startNode.
     * @return 
     */
    protected long createThe3NodeBitstring(int[] s) {
        assert(s.length == 3);
        long path = concatenate(0, 0, s);
        return path;
    }
    
    /**
     * for the pathNodeNumber which is the order number of the node in the path 
     * (e.g. second node in the path), set the node number to be the base10Node.
     * @param base10Node
     * @param path
     * @param pathNodeNumber number of the node within the path.  NOTE: this excludes
     * start node 0 and end node 0 (pathNodeNumber=0 corresponds to the
     * 2nd node in the final solution for the completed path for the algorithm instance).
     * @return 
     */
    protected long setBits(final int base10Node, final long path, final int pathNodeNumber) {
        long bitstring = path;
        long b = pathNodeNumber*w;
        long end = b + w;
        int b0;
        for (b = pathNodeNumber*w, b0=0; b < end; ++b, b0++) {
            if ((base10Node & (1L << b0)) != 0) { // test bit b0 in base10Node
                bitstring |= (1L << b); // set bit b in bitstring
            }
        }
        return bitstring;
    }
    
    /**
     * @param path encoded bit-string of ordered path nodes used in the memo.
     * @param nPathNodesSet the number of nodes currently set in the path
     * @param base10Nodes base 10 node indexes in the order to set into the bit-string path.
     * NOTE that s should not contain the startNode.
     */
    protected long concatenate(long path, int nPathNodesSet, int[] base10Nodes) {
        assert(numberOfSetNodes(path) == nPathNodesSet);
        long path2 = path;
        int i, si;
        for (i = 0; i < base10Nodes.length; ++i) {
            si = base10Nodes[i];
            assert(si != startNode);
            path2 = setBits(si, path2, nPathNodesSet + i);
        }
        return path2;
    }
    
    protected void compareToMin(long path, long sum) {
        assert(numberOfSetNodes(path) == (dist.length-1));
        int node1 = getBase10NodeIndex(0, path);
        int noden1 = getBase10NodeIndex(dist.length - 2, path);

        int ends = dist[startNode][node1] + dist[noden1][startNode];
        
        sum += ends;
        
        if (sum == minCost) {
            minPath.add(path);
        } else if (sum < minCost) {
            minCost = sum;
            minPath.clear();
            minPath.add(path);
        }
    }
    
    public long getMinCost() {
        return this.minCost;
    }
    public TLongList getMinPaths() {
        return new TLongArrayList(minPath);
    }

    /**
     * 
     * @param path a memo bit-string
     * @return 
     */
    protected int numberOfUnsetNodes(long path) {
        // composed of w-bit bit-strings
        int nn = (dist.length - 1) * w;
        
        // the first 64-((dist.length-1*w) can be discarded as unused space
        int nd = 64 - nn;
        int n0 = Long.numberOfLeadingZeros(path);
        assert(n0 >= nd);
        int nc = n0 - nd;
        int nUnset = nc/w;
        return nUnset;
    }
    
    /**
     * 
     * @param path a memo bit-string
     * @return 
     */
    protected int numberOfSetNodes(long path) {
        int nUnset = numberOfUnsetNodes(path);
        int nSet = (dist.length - 1) - nUnset;
        return nSet;
    }
}
