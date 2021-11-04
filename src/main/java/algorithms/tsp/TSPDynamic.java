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
 <pre>
 A completely dynamic solution requires a large amount of memory (see statiement below).
 Here is an outline of one, based upon what I learned from making the hybrid
 dynamic and brute force class:
 
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

The


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
public class TSPDynamic extends AbstractTSP {
    
    public TSPDynamic(double[][] dist) {
        super(dist);
    }
}
