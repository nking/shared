package algorithms.tsp;

import algorithms.Permutations;
import algorithms.PermutationsWithAwait;
import algorithms.SubsetChooser;
import algorithms.misc.MiscMath0;
import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Stack;

/**
 <pre>
 A completely dynamic solution requires a large amount of memory (see statiement below).
 Here is an outline of one, based upon what I learned from making the hybrid
 dynamic and brute force class:
 
 example for k=3, though will use k=2 in implementation:
 
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
need (n-1) > k for each stage.

  so consider (n-1) > k_prev and (n-1) .leq. k_next, cannot complete the next full
     permutation, so will finish the end of the paths using subsequences of size 3
     and without storing the partial results as they would mostly not be re-usable
     for the decreasing path lengths.  this separate section of iteration is
     like the brute-force section of TSPHybridDynamicBruteForce.java

  e.g.  consider (n-1) > 12 and (n-1) .leq. 24
                  n > 13 and n .leq. 25
  the permutations in between k=13 and k=25 can be handled in blocks of descending k.

  then use iterations over k-=3 of subsets and permutations, w/ a stack to
     finish complete paths and compare the completed path to minCost.
     
For the dynamic approach just outlined in the (n-1) > k section where k=3:
    n    c             compare to n!
    8    210           40320
   14    1.24e6        8.7e10
   29    1.46e16       8.8e31
   49    3.3e20 
   50    3.9e38
  731    5.58e538      7.7e1778
 where n is the number of nodes (=dist.length),
 and c is the number of elements one needs to store in a memo datastructure.

 For  (n-1) > k  where k=2:
    n    c             
    8    42           
   14    1.7e4       
   29    1.3e12       
   49    4.7e27 
   50    3.9e38
  731    6.6e711       
</pre>

see AbstractTSP for details about the keys used for memo.

 * @author nichole
 */
public class TSPDynamic extends AbstractTSP {
    
    public TSPDynamic(double[][] dist) {
        super(dist);
    }
    
    public void solve() throws InterruptedException {
        
        if (true) {
            throw new UnsupportedOperationException("not yet implemented");
        }
        
        if (minCost != sentinel) {
            reset();
        }
        int k0 = 2;
        int k = k0;
        int n = dist.length;
        
        int nNodesSet = k;
                 
        if ((dist.length <= (k + 1))) {
            initNodePaths();
            nNodesSet = dist.length - 1;
        } else {
            //dist.length > (k+1)
            if (k == 2) {
                // if k=2, init the k=4 paths in memo
                init4NodePaths();
            } else {
                // k == 3, init the k=3 paths in memo
                init3NodePaths();
            }
        }
        
        int nNodesRemaining = (n-1) - nNodesSet;
        
        // for the k increasing paths
        Stack<StackP> stack = new Stack<StackP>();
        
        // for the k too large paths
        Stack<StackP> stackDecr = new Stack<StackP>();
        
        // visit the initial path nodes in memo
        TLongDoubleIterator iter = memo.iterator();
        long bitstring, bitstring2;
        double sum, sum2;
        int nNodesRemaining2, k2;
        StackP currentStackP;
                
        for (int i = 0; i < memo.size(); ++i) {
            iter.advance();
            bitstring = iter.key();
            sum = iter.value();

            assert (stack.isEmpty());

            stack.add(new StackP(bitstring, sum, nNodesRemaining));

            while (!stack.isEmpty()) {

                currentStackP = stack.pop();
                bitstring2 = currentStackP.bitstring;
                sum2 = currentStackP.sum;
                nNodesRemaining2 = currentStackP.nNodesRemaining;
                k2 = dist.length - 1 - nNodesRemaining2;// number of bits in path bitstring2
                assert(k2 == this.numberOfSetNodes(bitstring2));
                
                if (nNodesRemaining2 == 0) {
                    compareToMin(bitstring2, sum2);
                    continue;
                }
                
                TIntList remaining = new TIntArrayList();
                findUnsetBitsBase10(bitstring2, remaining);
                assert (nNodesRemaining2 == remaining.size());
                int nBitsSet2 = (dist.length - 1) - nNodesRemaining2;
                                
                if (nNodesRemaining2 <= k2) {
                    stackDecr.add(new StackP(bitstring2, sum2, nNodesRemaining2));                 
                    continue;
                }
                                
                int lastNode = getBase10NodeIndex(nBitsSet2-1, bitstring2);

                // there are more than 3 nodes not yet set so can use subsetchooser
                SubsetChooser chooser = new SubsetChooser(nNodesRemaining2, k2);
                int[] sel = new int[k2];
                int s3, i3;
                int np = (int)MiscMath0.factorial(k2);
                
                int[] selPerm = new int[k2];
                double sum3;
                long path3, perm3i;
                int[] sel3 = new int[k2];
                PermutationsWithAwait permutations;
                while (true) {
                    s3 = chooser.getNextSubset(sel);
                    if (s3 == -1) {
                        break;
                    }

                    //transform sel to the bitstring2 unset indexes
                    for (i3 = 0; i3 < k2; ++i3) {
                        sel3[i3] = remaining.get(sel[i3]);
                    }
                    
                    permutations = new PermutationsWithAwait(sel3);

                    for (i3 = 0; i3 < np; ++i3) {
                        
                        permutations.getNext(selPerm);
                        
                        perm3i = createAMemoNodeBitstring(selPerm);
                        assert(memo.containsKey(perm3i));

                        // the number of set bits in path perm3i is the same as 
                        //    the number of set bits in bitstring2
                        //    and both are in memo already.
                        // path3 = bitstring2 + edge + perm31                                
                        sum3 = sum2 + dist[lastNode][selPerm[0]] + memo.get(perm3i);

                        path3 = concatenate(bitstring2, nBitsSet2, selPerm);

                        memo.put(path3, sum3);
                        
                        stack.add(new StackP(path3, sum3, nNodesRemaining2 - k2));
                    }
                }
            }
        }
        
        while (!stackDecr.isEmpty()) {
            //impl
        }
    }
}
