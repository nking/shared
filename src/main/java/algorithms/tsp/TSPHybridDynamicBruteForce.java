package algorithms.tsp;

import algorithms.Permutations;
import algorithms.SubsetChooser;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Stack;

/**
 Keeping notes here on looking at bit-string patterns to find ways to 
 compress the sums of 3-paths used in the subproblems
 for a dynamic solution to TSP as an exercise in recursion and
 dynamic programming (not a preferred impl due to the large runtime complexity). 

<pre>
    let start node = 0.  n=5
    The total number of permutations of a path with fixed start node = 0 and n=5
    is n!/n = 24.

    for n=127, nPerm = 2.4e+211, so if wanted to store each permutation, indexed,
    would need to use java's BigInteger as a bitstring or this project's VeryLongBitString.java
    (Not pursuing this or hash function options at this time).

    To re-use solved sub-problems, can determine how many sets of paths of 3 nodes are in
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
    (can use SubsetChooser.java and one of the 2 Permutation*.java classes).)

    n=10, k=3 .geq. count = 84

    Once all of the 3-node path subsets are generated and stored,
    one can consider again that each of the complete permutations of paths
    is composed of disjoint combinations of (n-1)/3 of the 3-node path subsets,
    (plus up to 2 nodes.  e.g. if n = 12, we have (n-1)/3 3-node path subsets
     + 1 fixed node + 2 free nodes that belong in the permuation).

    This code uses memoization for all 3-node combinations, then for each of those in
    memo, it uses brute force for the remaining path building where the brute force
    approach continues to use subsets and permutations and cost sums to complete
    each path and compare it to the min cost without further storing partial results in memo..

    see AbstractTSP class comments for notes about the keys used for memo.
          
    This class contains a recursive and non-recursive version of the same hybrid dynamic and brute
    force algorithm.
           
 * </pre>
 * 
 * @author nichole
 */
public class TSPHybridDynamicBruteForce extends AbstractTSP {
   
    /**
     *
     @param dist
     */
    public TSPHybridDynamicBruteForce(double[][] dist) {
        super(dist);
    }
    
    /**
     * this version is still roughly factorial.  its re-use of solving sub-problems
     * is only for the first 3 nodes in each path.
     * @throws java.lang.InterruptedException
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
        long[] memoKeys = memo.keys();
        long bitstring;
        double sum;
        for (int i = 0; i < memo.size(); ++i) {
            bitstring = memoKeys[i];
            sum = memo.get(bitstring);
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
     * is only for the first 3 nodes in each path.
     * @throws java.lang.InterruptedException
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
        long[] memoKeys = memo.keys();
        long bitstring, bitstring2;
        double sum, sum2;
        int nNodesRemaining2;
        StackP currentStackP;
        
        boolean storeInMemo = false;
        
        for (int i = 0; i < memo.size(); ++i) {
            bitstring = memoKeys[i];
            sum = memo.get(bitstring);

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
                                
                createAndStackSubsetPermutations(bitstring2, sum2, nBitsSet2, k, stack, storeInMemo);                
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
       
    /**
     *
     @return
     */
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
