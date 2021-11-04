package algorithms.tsp;

import algorithms.Permutations;
import algorithms.SubsetChooser;
import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Stack;

/**
 Keeping notes here on looking at bit-string patterns to find ways to 
 condense the sums of 3-paths used in the subproblems
 for a dynamic solution to TSP as an exercise in recursion and
 dynamic programming (not a preferred impl due to the large runtime complexity). 

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

    see AbstractTSP class comments for notes about the keys used for memo.

          
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
public class TSPHybridDynamicBruteForce extends AbstractTSP {
   
    
    public TSPHybridDynamicBruteForce(double[][] dist) {
        super(dist);
    }
    
    /**
     * this version is still roughly factorial.  its re-use of solving sub-problems
     * is only for the first 3 nodes.
     */
    public void solveRecursively() throws InterruptedException {
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
   
    /**
     * this version is still roughly factorial.  its re-use of solving sub-problems
     * is only for the first 3 nodes.
     */
    public void solveIteratively() throws InterruptedException {
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
                
                int nBitsSet2 = (dist.length - 1) - nNodesRemaining2;
                
                assert (nNodesRemaining2 == remaining.size());
                               
                if (nNodesRemaining2 <= k) {
                    
                    createAndStackPermutations(bitstring2, sum2,
                        nBitsSet2, remaining, stack, false);
                    
                    continue;
                }
                   
                
                int lastNode = getBase10NodeIndex(nBitsSet2-1, bitstring2);
          
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
                int[] sel3 = new int[k];
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
                        perm3i = createAMemoNodeBitstring(selPerm[i3]);
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
       
    public TLongList getMinPathBitstrings() {
        return new TLongArrayList(minPath);
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
                
                perm3i = createAMemoNodeBitstring(selPerm[i]);
                assert(memo.containsKey(perm3i));
                
                sum2 = sum + dist[lastNode][selPerm[i][0]] + memo.get(perm3i);
                
                path2 = concatenate(bitstring, nBitsSet, selPerm[i]);

                // no need to store                
                //memo.put(path2, sum2);
                
                r3(path2, sum2, nNodesRemaining - k);
            }
        }  
    }        

}
