package algorithms.tsp;

import algorithms.Permutations;
import algorithms.SubsetChooser;
import algorithms.misc.MiscMath0;
import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.TLongDoubleMap;
import gnu.trove.map.hash.TLongDoubleHashMap;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Stack;

/**
 * * <pre>
 * 
 * Keeping notes here on looking at bitstring patterns to find ways to 
 * condense the sums of 3-paths used in the subproblems
 * for a dynamic solution to TSP as an exercise in recursion and
 * dynamic programming (not a preferred impl due to the large runtime complexity).
 * 
 <pre>
 A completely dynamic solution requires a large amount of memory (see statiement below).
 Here is an outline of one, based upon what I learned from making the hybrid
 dynamic and brute force approach below:
 
 subset chooser for n, k w/ k=3
   permutations for k!
   store each k=3 permutation in memo
 subset chooser for n, k w/ k=6
   permutations for k!
   path should be composed of complete subparts + a connecting edge:
       p = p3_0 + p3_1 + edge where p3_0 and p3_1 are in memo from k=3 permutations
   store each k=6 permutation in memo
 subset chooser for n, k w/ k=12
   permutations for k!
   path should be composed of complete subparts + a connecting edge:
       p = p6_0 + p6_1 + edge where p6_0 and p6_1 are in memo from k=6 permutations
   store each k=12 permutation in memo
...
For the dynamic approach just outlined:
    n    c             compare to n!
    8    210           40320
   14    1.24e6        8.7e10
   29    1.46e16       8.8e31
  731    5.58e538      7.7e1778
 where n is the number of cities,
 and c is the number of elements one needs to store in a memo datastructure.
 
 </pre>
 * 
<pre>
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
               if (nNodesRemaining .leq. k) {
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
                  
                   then make the 6 mermutations of the 3 set bits
                   and for each of those:
                   
                       bitstring2 = concatenate(bitstring, nBitsSet, si);
                       sum2 = sum + memo.get(si);  then use a method to convert si to a memo bitstring
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
                   if (nNodesRemaining2 .leq. k) {
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
                       
                       then make the 6 mermutations of the 3 set bits
                       and for each of those:
                   
                           bitstring3 = concatenate(bitstring2, nBitsSet, si);
                           sum3 = sum2 + memo.get(si); then use method to convert si to a memo bitstring
                           memo.set(bitstring3, sum3);

                           stack.add(bitstring3, sum3, nNorderRemaining - 3);
                   }
               }
           }
    
 * </pre>
 * 
 * @author nichole
 */
public class TSPHybridDynamicBruteForce {
   
    private double minCost = Double.POSITIVE_INFINITY;
    private final TLongList minPath = new TLongArrayList();
    private final int startNode = 0;
    private final double[][] dist;
    //TODO: consider other formats for memo
    private final TLongDoubleMap memo;
    private final double sentinel = Double.POSITIVE_INFINITY;
    private final long totalNPerm;
    private final long totalNSubSet;
    private final long totalNSubSeq;
    private final int w; // number of bits a city takes in a path where path is a bitstring of type long
    
    public TSPHybridDynamicBruteForce(double[][] dist) {
        
        this.dist = dist;
        int n = dist.length;
        long nPerm = MiscMath0.factorial(n); // max for n=13 for limit of array length
        totalNPerm = nPerm/n;
        
        //TODO: add in the number of permutations for those not in a 3-set, that is,
        //   the 2 node and 1-node permutations
        totalNSubSet = countTotalNumSubSetInvocations(n - 1); // max for n=338 for limit of array length
        totalNSubSeq = countTotalNumSubSeqInvocations(n - 1); 
        
        System.out.printf("nPerm=%d, totalNSubSet=%d  totalNSubSeq=%d\n", 
            totalNPerm, totalNSubSet, totalNSubSeq);
        
        int sz = (int)MiscMath0.computeNDivNMinusK(dist.length-1, 3);
        
        /*if (totalNSubSeq > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("this class can solve for 13 cities at most."
                + " one could design a datastructure to hold more than (1<<31)-1 items, "
                + " but the algorithms in this class have exponential runtime complexity,"
                + " so it would be better to choose a good approximate TSP.");
        }*/
        n = (int)totalNPerm;
                        
        w = (int)(Math.ceil(Math.log(dist.length)/Math.log(2)));
        
        memo = new TLongDoubleHashMap(sz);
        
        /*
        //find max n cities for limit of array length
        int nc = 100;
        long wc, nb;
        while (true) {
            long c = countTotalNumSubSeqInvocations(nc);
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
    
    private void reset() {
        minCost = sentinel;
        minPath.clear();
        memo.clear();
    }
    
    /**
     * this version is still roughly factorial.  its re-use of solving sub-problems
     * is only for the first 3 nodes.
     */
    public void solveRecursively() {
        if (minCost != sentinel) {
            reset();
        }
        
        int n = dist.length;
        int nNodesSet = 3;
               
        // initialize memo with the first 3-node paths to re-use in permuations
        
        if (dist.length <= (3 + 1)) {
            // cannot use subsetchooser, so go straight to permutations
            initNodePaths();
            nNodesSet = dist.length - 1;
        } else {
            init3NodePaths();
        }
        
        cr0 = memo.size();
        
        int nNodesRemaining = (n-1) - nNodesSet;
        
        // visit the initial path nodes in memo
        TLongDoubleIterator iter = memo.iterator();
        long bitstring;
        double sum;
        for (int i = 0; i < memo.size(); ++i) {
            iter.advance();
            bitstring = iter.key();
            sum = iter.value();
            r3(bitstring, sum, nNodesRemaining);
        }
        //r0=240 n=6
        //cr0=72 n=5
        //cr0=12 n=4
        double dyn = Math.pow(2, dist.length) + Math.pow(dist.length, 2);
        double dyn1 = Math.pow(2, dist.length-1) + Math.pow(dist.length-1, 2);
        System.out.printf("cr0=%d n=%d\n", cr0, dist.length);
        System.out.printf("    totalNPerm=%d totalNSubSet=%d totalNSubSeq=%d  dyn=%.1f dyn1=%.1f\n", 
            totalNPerm, totalNSubSet, totalNSubSeq, dyn, dyn1);
    }
    long cr0;
    
    private static class StackP {
        long bitstring;
        double sum;
        int nNodesRemaining;
        public StackP(long path, double cost, int nRemaining) {
            this.bitstring = path;
            this.sum = cost;
            this.nNodesRemaining = nRemaining;
        }
    }
    
    /**
     * this version is still roughly factorial.  its re-use of solving sub-problems
     * is only for the first 3 nodes.
     */
    public void solveIteratively() {
        if (minCost != sentinel) {
            reset();
        }
        
        int n = dist.length;
        int nNodesSet = 3;
               
        // initialize memo with the first 3-node paths to re-use in permuations
                
        if (dist.length <= (3 + 1)) {
            // cannot use subsetchooser, so go straight to permutations
            initNodePaths();
            nNodesSet = dist.length - 1;
        } else {
            init3NodePaths();
        }
        
        long c0 = memo.size();
        long c1 = memo.size();
        long c2 = memo.size();
  
        int nNodesRemaining = (n-1) - nNodesSet;
        
        Stack<StackP> stack = new Stack<StackP>();
        
        // visit the initial path nodes in memo
        TLongDoubleIterator iter = memo.iterator();
        long bitstring, bitstring2;
        double sum, sum2;
        int nNodesRemaining2;
        StackP currentStackP;
        
        for (int i = 0; i < memo.size(); ++i) {
            iter.advance();
            bitstring = iter.key();
            sum = iter.value();

            assert (stack.isEmpty());

            stack.add(new StackP(bitstring, sum, nNodesRemaining));

            while (!stack.isEmpty()) {
                c0++;

                currentStackP = stack.pop();
                bitstring2 = currentStackP.bitstring;
                sum2 = currentStackP.sum;
                nNodesRemaining2 = currentStackP.nNodesRemaining;
                
                if (nNodesRemaining2 == 0) {
                    compareToMin(bitstring2, sum2);
                    continue;
                }
                
                c1++;
                
                int k = 3;
                
                TIntList remaining = new TIntArrayList();
                findUnsetBitsBase10(bitstring2, remaining);
                assert (nNodesRemaining2 == remaining.size());
                int nBitsSet2 = (dist.length - 1) - nNodesRemaining2;
                
                int lastNode = getBase10NodeIndex(nBitsSet2-1, bitstring2);
                
                if (nNodesRemaining2 <= k) {
                    int firstNode = getBase10NodeIndex(0, bitstring2);

                    if (nNodesRemaining2 == 2) {
                        // 2 permutations, add each to the end of the path.
                        long bitstring3 = setBits(remaining.get(0), bitstring2, dist.length - 1 - 2);
                        bitstring3 = setBits(remaining.get(1), bitstring3, dist.length - 1 - 1);

                        long bitstring4 = setBits(remaining.get(1), bitstring2, dist.length - 1 - 2);
                        bitstring4 = setBits(remaining.get(0), bitstring4, dist.length - 1 - 1);

                        double sum3 = sum + dist[lastNode][remaining.get(0)] + 
                            dist[remaining.get(0)][remaining.get(1)];
                        double sum4 = sum + dist[lastNode][remaining.get(1)] + 
                            dist[remaining.get(1)][remaining.get(0)];

                        stack.add(new StackP(bitstring3, sum3, 0));
                        stack.add(new StackP(bitstring4, sum4, 0));
                    } else {
                        // 1 permutation, meaning, add the node to the end of the path
                        int node = remaining.get(0);
                        long bitstring3 = setBits(remaining.get(0),
                            bitstring2, dist.length - 1 - 1);
                        double sum3 = sum2 + dist[lastNode][remaining.get(0)];

                        // no need to save these
                        //memo.put(bitstring1, sum1);
                        stack.add(new StackP(bitstring3, sum3, 0));
                    }
                    continue;
                }
                
                // there are more than 3 nodes not yet set so can use subsetchooser
                SubsetChooser chooser = new SubsetChooser(nNodesRemaining2, k);
                int[] sel = new int[k];
                int s3, i3;
                int[][] selPerm = new int[6][k];
                for (i3 = 0; i3 < selPerm.length; ++i3) {
                    selPerm[i3] = new int[k];
                }
                double sum3;
                long path3, perm3i;
                int[] sel3 = new int[nNodesRemaining];
                while (true) {
                    s3 = chooser.getNextSubset(sel);
                    if (s3 == -1) {
                        break;
                    }

                    //transform sel to the bitstring2 unset indexes
                    for (i3 = 0; i3 < k; ++i3) {
                        sel3[i3] = remaining.get(sel[i3]);
                    }

                    Permutations.permute(sel3, selPerm);

                    for (i3 = 0; i3 < selPerm.length; ++i3) {
                        c2++;
                        perm3i = createThe3NodeBitstring(selPerm[i3]);
                        assert(memo.containsKey(perm3i));

                        sum3 = sum2 + dist[lastNode][selPerm[i3][0]] + memo.get(perm3i);

                        path3 = concatenate(bitstring2, nBitsSet2, selPerm[i3]);

                        // no need to store
                        //memo.put(path2, sum2);
                        stack.add(new StackP(path3, sum3, nNodesRemaining2 - k));
                    }
                }
            }
        } 
        //c0=240 c1=120 c2=60 n=6
        //c0=72 c1=48 c2=24 n=5
        //c0=12 c1=6 c2=6 n=4
        // if dynamic: O(n^2 * 2^n) 
        double dyn = Math.pow(2, dist.length) + Math.pow(dist.length, 2);
        double dyn1 = Math.pow(2, dist.length-1) + Math.pow(dist.length-1, 2);
        System.out.printf("c0=%d c1=%d c2=%d n=%d\n", c0, c1, c2, dist.length);
        System.out.printf("    totalNPerm=%d totalNSubSet=%d totalNSubSeq=%d  dyn=%.1f dyn1=%.1f\n", 
            totalNPerm, totalNSubSet, totalNSubSeq, dyn, dyn1);
    }
   
    // total number of subsetchooser invocations.  max n = 507 for count < Integer.MAX_VALUE
    protected long countTotalNumSubSeqInvocations(int n) {
        int k=3, n2=n;
        long c = 0;
        while (n2 > k) {
            c += MiscMath0.computeNDivKTimesNMinusK(n2, k);
            n2 -= k;
        }
        return c;
    }
    
    // total number of k-permutations, that is, n!/(n-k)!.  using the subset chooser
    //   to find the sets of size 3, then permuting each of those into
    //   distinct ordered sequences
    protected long countTotalNumSubSetInvocations(int n) {
        int k=3, n2=n;
        long c = 0;
        while (n2 > k) {
            c += MiscMath0.computeNDivNMinusK(n2, k);
            n2 -= k;
        }
        return c;
    }
    
    /**
     * roughly counting k-permutations for a dynamic approach where k is
     * increased by a factor of 2 each time and begins with k=3.
     * @param n
     * @return 
     */
    protected static BigInteger count0(int n) {
        n = n - 1;
        BigInteger c1, c2, c3;
        BigInteger c0 = BigInteger.ZERO;
        int k = 3;
        while ((n-k) > k) {
            c1 = MiscMath0.computeNDivKTimesNMinusKBigInteger(n, k);
            c2 = MiscMath0.factorialBigInteger(k);
            c3 = c1.multiply(c2);
            c0 = c0.add(c3);
            k *= 2;
        }
        return c0;
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
        long bitstring = path; // 11 10 01
        long b = pathNodeNumber*w;
        long end = b + w;
        int b0;
        for (b = pathNodeNumber*w, b0=0; b < end; ++b, b0++) {
            if ((base10Node & (1L << b0)) != 0) { // test bit b0 in pathNodeNumber
                bitstring |= (1L << b); // set bit b in bitstring
            } else {
                bitstring &= ~(1L << b);//clear bit
            }
        }
        return bitstring;
    }
    
    /**
     * read the given bit-string encoded for use with memo, to find the
     * set bits and return the nodes as a base10 bit-string (without the path order information).
     * @param bitstring
     * @return 
     */
    protected long findSetBitsBase10(long bitstring) {
        long base10nodes = 0;
        int i, b, node, j;
        int bf = w*(dist.length - 1);
        for (b = 0; b < bf; b+= w) {
            node = 0;
            for (i = b, j=0; i < (b + w); ++i, ++j) {
                if ((bitstring & (1L << i)) != 0) {
                    node += (1 << j);
                }
            }
            if (node > 0) {
                base10nodes |= (1L << node);
            }
        }
        return base10nodes;
    }
    
    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @return 
     */
    protected long findUnsetBitsBase10(long bitstring) {
        long base10nodesSet = findSetBitsBase10(bitstring);
        long base10NotSet = 0;
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) == 0) {
                base10NotSet |= (1L << b);
            }
        }
        return base10NotSet;
    }
    
    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @param out list to hold the output node indexes 
     */
    protected void findUnsetBitsBase10(long bitstring, TIntList out) {
        out.clear();
        long base10nodesSet = findSetBitsBase10(bitstring);
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) == 0) {
                out.add(b);
            }
        }
    }
    
    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string that has lost
     * information about the path node order.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @param out list to hold the output node indexes 
     */
    protected void findSetBitsBase10(long bitstring, TIntList out) {
        out.clear();
        long base10nodesSet = findSetBitsBase10(bitstring);
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) != 0) {
                out.add(b);
            }
        }
    }
    
    /**
     * read path into base10 node numbers, preserving order of path
     * @param bitstring
     * @param out 
     */
    protected void readPathIntoBase10(long bitstring, TIntList out) {
        out.clear();
        
        int i, b, node, j;
        int bf = w*(dist.length - 1);
        for (b = 0; b < bf; b+= w) {
            node = 0;
            for (i = b, j=0; i < (b + w); ++i, ++j) {
                if ((bitstring & (1L << i)) != 0) {
                    node += (1 << j);
                }
            }
            if (node > 0) {
                out.add(node);
            }
        }
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
    
    protected void compareToMin(long path, double sum) {
        assert(numberOfSetNodes(path) == (dist.length-1));
        
        int node1 = getBase10NodeIndex(0, path);
        int noden1 = getBase10NodeIndex(dist.length - 2, path);

        double ends = dist[startNode][node1] + dist[noden1][startNode];
        
        double sum2 = sum + ends;
        
        //debug
        /*TIntList p = new TIntArrayList();
        readPathIntoBase10(path, p);
        System.out.printf("final: bs=%s (%s) sum=%.2f sum2=%.2f, min=%.2f\n",
            Long.toBinaryString(path),
            Arrays.toString(p.toArray()), sum, sum2, minCost);
        */
        // end debug
        
        if (sum2 == minCost) {
            minPath.add(path);
        } else if (sum2 < minCost) {
            minCost = sum2;//3, 2, 4, 1, 5
            minPath.clear();
            minPath.add(path);
        }
    }
    
    public double getMinCost() {
        return this.minCost;
    }
    
    public int getNumberOfMinPaths() {
        return minPath.size();
    }
    
    public TIntList getMinPath(int idx) {
        long path = minPath.get(idx);
        TIntList out = new TIntArrayList();
        readPathIntoBase10(path, out);
        out.insert(0, startNode);
        out.add(startNode);
        return out;
    }
    
    public TLongList getMinPathBitstrings() {
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
    
    /**
     * initialize memo with permutations for all unset path nodes, where the
     * number of unset path nodes < 4.
     */
    protected void initNodePaths() {
        assert(memo.isEmpty());
        
        int nUnset = dist.length - 1;
        
        int i;
        int[] sel = new int[nUnset];
        for (i = 1; i <= sel.length; ++i) {
            sel[i-1] = i;
        }
        
        int nPerm = (int)MiscMath0.factorial(nUnset);
        
        final int[][] selPerm = new int[nPerm][nUnset];
        for (i = 0; i < selPerm.length; ++i) {
            selPerm[i] = new int[nUnset];
        }
        
        Permutations.permute(sel, selPerm);
        
        double sum;
        long path;
        int j, i0, i1;
        for (i = 0; i < selPerm.length; ++i) {
            sum = 0;
            
            //System.out.println("    selPerm=" + Arrays.toString(selPerm[i]));

            path = createThe3NodeBitstring(selPerm[i]);

            for (j = 1; j < selPerm[i].length; ++j) {
                i0 = selPerm[i][j - 1];
                i1 = selPerm[i][j];
                sum += dist[i0][i1];
            }
            memo.put(path, sum);
        }
    }

    /**
     * initialize memo with permutations for all subsets of 3 path nodes, where the
     * number of unset path nodes is > 3.
     */
    protected void init3NodePaths() {        
        assert(memo.isEmpty());
        
        int k = 3;
        final int[] sel = new int[k];
        final int[] sel2 = new int[k];
        int s, i;
        final int[][] selPerm = new int[6][k];
        for (i = 0; i < selPerm.length; ++i) {
            selPerm[i] = new int[k];
        }
        
        TIntList remaining = new TIntArrayList();
        findUnsetBitsBase10(0, remaining);
        //System.out.println("remaining unset=" + remaining.toString());

        int j, i0, i1;
        long path, sum;
        SubsetChooser chooser = new SubsetChooser(dist.length-1, k);
        while (true) {
            s = chooser.getNextSubset(sel);
            if (s == -1) {
                break;
            }
            //System.out.println("sel=" + Arrays.toString(sel));
            //transform sel to the bitstring unset indexes
            for (i = 0; i < k; ++i) {
                sel2[i] = remaining.get(sel[i]);
            }
            //System.out.println("    sel2=" + Arrays.toString(sel2));

            Permutations.permute(sel2, selPerm);
            
            for (i = 0; i < selPerm.length; ++i) {
                sum = 0;
                path = createThe3NodeBitstring(selPerm[i]);
                
                for (j = 1; j < k; ++j) {
                    i0 = selPerm[i][j-1];
                    i1 = selPerm[i][j];
                    sum += dist[i0][i1];
                }
                memo.put(path, sum);
            }
        }
    }
    
    private void r3(long bitstring, double sum, int nNodesRemaining) {
        cr0++;
        //debug
        /*
        TIntList p = new TIntArrayList();
        readPathIntoBase10(bitstring, p);
        System.out.printf("bs=%s (%s) sum=%.2f nR=%d\n",
            Long.toBinaryString(bitstring),
            Arrays.toString(p.toArray()), sum, nNodesRemaining);
        */
        // end debug
        
        if (nNodesRemaining == 0) {
            compareToMin(bitstring, sum);
            return;
        }
        
        int k = 3;
        
        TIntList remaining = new TIntArrayList();
        findUnsetBitsBase10(bitstring, remaining);
        assert (nNodesRemaining == remaining.size());
        int nBitsSet = (dist.length - 1) - nNodesRemaining;
        
        if (nNodesRemaining <= k) {
            int firstNode = getBase10NodeIndex(0, bitstring);
            int lastNode = getBase10NodeIndex(nBitsSet-1, bitstring);
            
            if (nNodesRemaining == 2) {
                
                // 2 permutations, add each to the end of the path.
                long bitstring1 = setBits(remaining.get(0), bitstring, dist.length - 1 - 2);
                bitstring1 = setBits(remaining.get(1), bitstring1, dist.length - 1 - 1);

                long bitstring2 = setBits(remaining.get(1), bitstring, dist.length - 1 - 2);
                bitstring2 = setBits(remaining.get(0), bitstring2, dist.length - 1 - 1);

                double sum1 = sum + dist[lastNode][remaining.get(0)] + 
                    dist[remaining.get(0)][remaining.get(1)];
                double sum2 = sum + dist[lastNode][remaining.get(1)] + 
                    dist[remaining.get(1)][remaining.get(0)];

                // no need to save these       
                //memo.put(bitstring1, sum1);
                //memo.put(bitstring2, sum2);

                r3(bitstring1, sum1, 0);
                r3(bitstring2, sum2, 0);
            } else {
                // 1 permutation, meaning, add the node to the end of the path
                int node = remaining.get(0);
                long bitstring1 = setBits(remaining.get(0), 
                    bitstring, dist.length - 1 - 1);            
                double sum1 = sum + dist[lastNode][remaining.get(0)];
                
                // debug
                int n1 = numberOfSetNodes(bitstring);
                int n2 = numberOfSetNodes(bitstring1);
                boolean t1 = (n1 == nBitsSet);
                long u2 = findUnsetBitsBase10(bitstring1);
                
                // no need to save these       
                //memo.put(bitstring1, sum1);
                r3(bitstring1, sum1, 0);
            }
            return;
        }
        
        int lastNode = getBase10NodeIndex(nBitsSet-1, bitstring);

        SubsetChooser chooser = new SubsetChooser(nNodesRemaining, k);
        int[] sel = new int[k];
        int s, i;
        int[][] selPerm = new int[6][k];
        for (i = 0; i < selPerm.length; ++i) {
            selPerm[i] = new int[k];
        }
        double sum2;
        long path2, perm3i;
        int[] sel2 = new int[nNodesRemaining];
        while (true) {
            s = chooser.getNextSubset(sel);
            if (s == -1) {
                break;
            }
            
            //transform sel to the bitstring unset indexes
            for (i = 0; i < k; ++i) {
                sel2[i] = remaining.get(sel[i]);
            }
            
            Permutations.permute(sel2, selPerm);
            
            for (i = 0; i < selPerm.length; ++i) {
                
                perm3i = createThe3NodeBitstring(selPerm[i]);
                assert(memo.containsKey(perm3i));
                
                sum2 = sum + dist[lastNode][selPerm[i][0]] + memo.get(perm3i);
                
                path2 = concatenate(bitstring, nBitsSet, selPerm[i]);

                // no need to store                
                //memo.put(path2, sum2);
                
                r3(path2, sum2, nNodesRemaining - k);
            }
        }  
    }        

    public int getMemoLength() {
        return memo.size();
    }

    protected void printMemo() {
        TLongDoubleIterator iter = memo.iterator();
        long bitstring;
        double sum;
        TIntList p = new TIntArrayList();
        
        for (int i = 0; i < memo.size(); ++i) {
            iter.advance();
            bitstring = iter.key();
            sum = iter.value();
            
            p.clear();
            readPathIntoBase10(bitstring, p);
            
            System.out.printf("memo: (%s) sum=%.2f min=%.2f\n",
                Arrays.toString(p.toArray()), sum, minCost);
        }
        System.out.flush();
    }
}
